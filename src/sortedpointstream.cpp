#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <memory>
#include <cstdio>
#include <cerrno>

#include "omp.h"

#include "liblas/liblas.hpp"

#include "geotools.h"
#include "sortedpointstream.hpp"
#include "lasutil.hpp"

using namespace geotools::las;
using namespace geotools::util;

double LASPoint::scaleX = 0;
double LASPoint::scaleZ = 0;
double LASPoint::scaleY = 0;

std::string _rowFile(uint32_t row) {
	std::stringstream fs;
	fs << "/tmp/sps_" << row << ".tmp";
	return fs.str();
}

uint32_t _getRow(double y, const Bounds &bounds, double blockSize) {
	if(blockSize < 0) {
		return (uint32_t) ((y - bounds.maxy()) / blockSize);
	} else {
		return (uint32_t) ((bounds.miny() - y) / blockSize);
	}
}


LASPoint::LASPoint() {}
LASPoint::LASPoint(const liblas::Point &pt) {
	update(pt);
}

void LASPoint::setScale(double x, double y, double z) {
	scaleX = x;
	scaleY = y;
	scaleZ = z;
}

void LASPoint::update(const liblas::Point &pt) {
	x = pt.GetX();
	y = pt.GetY();
	z = pt.GetZ();
	intensity = pt.GetIntensity();
	ret = pt.GetReturnNumber();
	numRets = pt.GetNumberOfReturns();
	cls = pt.GetClassification().GetClass();
	angle = pt.GetScanAngleRank();
}

void LASPoint::read(std::istream &str) {
	int32_t xx, yy, zz;
	str.read((char *) &xx, sizeof(int32_t));
	str.read((char *) &yy, sizeof(int32_t));
	str.read((char *) &zz, sizeof(int32_t));
	str.read((char *) &intensity, sizeof(uint16_t));
	str.read((char *) &ret, sizeof(uint16_t));
	str.read((char *) &numRets, sizeof(uint16_t));
	str.read((char *) &cls, sizeof(uint8_t));
	str.read((char *) &angle, sizeof(int8_t));
	x = (double) (xx * scaleX);
	y = (double) (yy * scaleY);
	z = (double) (zz * scaleZ);
}

void LASPoint::write(std::ostream &str) const {
	int32_t xx = (int32_t) (x / scaleX);
	int32_t yy = (int32_t) (y / scaleY);
	int32_t zz = (int32_t) (z / scaleZ);
	str.write((char *) &xx, sizeof(int32_t));
	str.write((char *) &yy, sizeof(int32_t));
	str.write((char *) &zz, sizeof(int32_t));
	str.write((char *) &intensity, sizeof(uint16_t));
	str.write((char *) &ret, sizeof(uint16_t));
	str.write((char *) &numRets, sizeof(uint16_t));
	str.write((char *) &cls, sizeof(uint8_t));
	str.write((char *) &angle, sizeof(int8_t));
}

void LASPoint::read(std::FILE *str) {
	int32_t xx, yy, zz;
	if(!std::fread((void *) &xx, sizeof(int32_t), 1, str))
		g_runerr("Failed to read x.")
	if(!std::fread((void *) &yy, sizeof(int32_t), 1, str))
		g_runerr("Failed to read x.")
	if(!std::fread((void *) &zz, sizeof(int32_t), 1, str))
		g_runerr("Failed to read x.")
	if(!std::fread((void *) &intensity, sizeof(uint16_t), 1, str))
		g_runerr("Failed to read x.")
	if(!std::fread((void *) &ret, sizeof(uint16_t), 1, str))
		g_runerr("Failed to read x.")
	if(!std::fread((void *) &numRets, sizeof(uint16_t), 1, str))
		g_runerr("Failed to read x.")
	if(!std::fread((void *) &cls, sizeof(uint8_t), 1, str))
		g_runerr("Failed to read x.")
	if(!std::fread((void *) &angle, sizeof(int8_t), 1, str))
		g_runerr("Failed to read x.")
	x = (double) (xx * scaleX);
	y = (double) (yy * scaleY);
	z = (double) (zz * scaleZ);
}

void LASPoint::write(std::FILE *str) const {
	int32_t xx = (int32_t) (x / scaleX);
	int32_t yy = (int32_t) (y / scaleY);
	int32_t zz = (int32_t) (z / scaleZ);
	std::fwrite((const void *) &xx, sizeof(int32_t), 1, str);
	std::fwrite((const void *) &yy, sizeof(int32_t), 1, str);
	std::fwrite((const void *) &zz, sizeof(int32_t), 1, str);
	std::fwrite((const void *) &intensity, sizeof(uint16_t), 1, str);
	std::fwrite((const void *) &ret, sizeof(uint16_t), 1, str);
	std::fwrite((const void *) &numRets, sizeof(uint16_t), 1, str);
	std::fwrite((const void *) &cls, sizeof(uint8_t), 1, str);
	std::fwrite((const void *) &angle, sizeof(int8_t), 1, str);
}

void LASPoint::read(void *str) {
	char *ptr = (char *) str;
	x         = *((int32_t *)  ptr) * scaleX; ptr += sizeof(int32_t);
	y         = *((int32_t *)  ptr) * scaleY; ptr += sizeof(int32_t);
	z         = *((int32_t *)  ptr) * scaleZ; ptr += sizeof(int32_t);
	intensity = *((uint16_t *) ptr);          ptr += sizeof(uint16_t);
	ret       = *((uint16_t *) ptr);          ptr += sizeof(uint16_t);
	numRets   = *((uint16_t *) ptr);          ptr += sizeof(uint16_t);
	cls       = *((uint8_t *)  ptr);          ptr += sizeof(uint8_t);
	angle     = *((int8_t *)   ptr);
}

void LASPoint::write(void *str) const {
	char *ptr = (char *) str;
	*((int32_t *)  ptr) = (int32_t) (x / scaleX); ptr += sizeof(int32_t);
	*((int32_t *)  ptr) = (int32_t) (y / scaleY); ptr += sizeof(int32_t);
	*((int32_t *)  ptr) = (int32_t) (z / scaleZ); ptr += sizeof(int32_t);
	*((uint16_t *) ptr) = intensity;              ptr += sizeof(uint16_t);
	*((uint16_t *) ptr) = ret;                    ptr += sizeof(uint16_t);
	*((uint16_t *) ptr) = numRets;                ptr += sizeof(uint16_t);
	*((uint8_t *)  ptr) = cls;                    ptr += sizeof(uint8_t);
	*((int8_t *)   ptr) = angle;
}

bool LASPoint::last() const {
	return numRets > 0 && ret == numRets;
}

bool LASPoint::first() const {
	return numRets > 0 && ret == 1;
}

bool LASPoint::intermediate() const {
	return numRets > 2 && ret > 1 && ret < numRets;
}

bool LASPoint::ground() const {
	return cls == 2;
}

bool LASPoint::single() const {
	return numRets == 1;
}


SortedPointStream::SortedPointStream(const std::list<std::string> &files, double blockSize, bool rebuild) :
	m_files(files),
	m_blockSize(blockSize),
	m_rebuild(rebuild) {
}

SortedPointStream::~SortedPointStream() {
	delete m_pf;
}

uint64_t SortedPointStream::pointCount() const {
	return m_pf->pointCount();
}

uint32_t SortedPointStream::rowCount() const {
	return m_pf->rowCount();
}

void SortedPointStream::init() {

	if(m_inited) return;
	m_inited = true;

	m_bounds.collapse();
	m_pointCount = 0;
	m_row = 0;

	liblas::ReaderFactory rf;
	std::vector<std::string> files0;
	Bounds workBounds; // Configure

	for(const std::string &file : m_files) {
		g_debug(" -- init: opening file: " << file);
		std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
		liblas::Reader lasReader = rf.CreateWithStream(instr);
		liblas::Header lasHeader = lasReader.GetHeader();
		LASPoint::setScale(lasHeader.GetScaleX(), lasHeader.GetScaleY(), lasHeader.GetScaleZ());
		g_debug(" -- init computing bounds");
		Bounds fileBounds;
		if(!LasUtil::computeLasBounds(lasHeader, fileBounds, 2))
			LasUtil::computeLasBounds(lasReader, fileBounds, 2); // If the header bounds are bogus.
		if(fileBounds.intersects(workBounds)) {
			m_bounds.extend(fileBounds);
			m_pointCount += lasHeader.GetPointRecordsCount();
			files0.push_back(file);
		}
	}
	m_bounds.snap(g_abs(m_blockSize));

	if(m_rebuild) {
		
		uint32_t rows = (uint32_t) (m_bounds.height() / m_blockSize) + 1;
		uint32_t rowLen = m_pointCount / rows * 2;
		uint32_t rowSize = (rowLen * LASPoint::dataSize) + (3 * sizeof(uint32_t));
		uint64_t size = rows * rowSize;
		uint32_t nextRow = 0;

		std::unique_ptr<MappedFile> mf = Util::mapFile(Util::tmpFile(), size);
		char* data = (char *) mf->data();

		//#pragma omp parallel for
		for(uint32_t i = 0; i< files0.size(); ++i) {
			
			const std::string &file = files0[i];
			g_debug(" -- init: opening file: " << file << "; " << omp_get_thread_num());
			
			std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
			liblas::Reader lasReader = rf.CreateWithStream(instr);
			liblas::Header lasHeader = lasReader.GetHeader();

			uint32_t row;
			uint32_t count;
			uint32_t jump;
			uint32_t offset;

			LASPoint pt;
			while(lasReader.ReadNextPoint()) {
				pt.update(lasReader.GetPoint());
				
				if(m_blockSize < 0) {
					row = (uint32_t) ((pt.y - m_bounds.maxy()) / m_blockSize);
				} else {
					row = (uint32_t) ((pt.y - m_bounds.miny()) / m_blockSize);
				}

				offset = row * rowSize + sizeof(uint32_t);
				if(offset >= size) {
					g_warn("Offset larger than size " << offset << ", " << size);
					break;
				}

				char *d = data + offset; // skip the row column


				jump = 0;
				do {
					offset = jump * rowSize + sizeof(uint32_t);
					if(offset >= size) {
						g_warn("Offset larger than size " << offset << ", " << size);
						break;
					}
					d = data + offset; // skip the row column
					count = *((uint32_t *) d); d += sizeof(uint32_t);
					jump = *((uint32_t *) d); d += sizeof(uint32_t);
				} while(jump > 0);

				if(count == rowLen) {
					++nextRow;
					jump = row + nextRow;
					offset = jump * rowSize + sizeof(uint32_t);
					if(offset >= size) {
						g_warn("Offset larger than size " << offset << ", " << size);
						break;
					}
					d -= sizeof(uint32_t) * 2; // back to the jump column
					*((uint32_t *) d) = jump;
					d = data + offset; // skip the row column
					jump = 0;
					count = 0;
				}


				*((uint32_t *) d) = row;       d += sizeof(uint32_t);
				*((uint32_t *) d) = count + 1; d += sizeof(uint32_t);
				*((uint32_t *) d) = 0;         d += sizeof(uint32_t);

				d += count * LASPoint::dataSize;

				++count;

				pt.write((void *) d);

			}
		}
	}
}

const Bounds& SortedPointStream::bounds() const {
	return m_bounds;
}

bool SortedPointStream::next(std::list<std::shared_ptr<LASPoint> > &pts) {
	return false;
}
