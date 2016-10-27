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
	m_rebuild(rebuild),
	m_inited(false),
	m_file(nullptr) {
}

SortedPointStream::~SortedPointStream() {
	if(m_file)
		std::fclose(m_file);
}

uint64_t SortedPointStream::pointCount() const {
	return m_pointCount;
}

uint32_t SortedPointStream::rowCount() const {
	return m_rowCount;
}

void _produce(SortedPointStream *s) {
	s->produce();
}

void SortedPointStream::produce() {

	liblas::ReaderFactory rf;
	while(m_fileq.size()) {

		std::string file;

		m_fmtx.lock();
		if(m_fileq.size()) {
			file = m_fileq.front();
			m_fileq.pop();
		}
		m_fmtx.unlock();

		std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
		liblas::Reader lasReader = rf.CreateWithStream(instr);
		liblas::Header lasHeader = lasReader.GetHeader();

		uint32_t row;
		LASPoint pt;
		while(lasReader.ReadNextPoint()) {
			pt.update(lasReader.GetPoint());
			
			if(m_blockSize < 0) {
				row = (uint32_t) ((pt.y - m_bounds.maxy()) / m_blockSize);
			} else {
				row = (uint32_t) ((pt.y - m_bounds.miny()) / m_blockSize);
			}

			m_cmtx.lock();
			m_cache[row].push_back(new LASPoint(pt));
			m_cmtx.unlock();
		}

	}
}

void _consume(SortedPointStream *s) {
	s->consume();
}

void SortedPointStream::consume() {

	uint32_t rowSize = 3 * sizeof(uint32_t) + 1024 * LASPoint::dataSize;

	while(m_running) {
		for(const auto &it : m_cache) {
			if(it.second.size() > 1024 || m_finalize) {

				uint32_t row = it.first;
				uint32_t count = g_min(1024, it.second.size());
				uint32_t curIdx = 0;
				uint32_t prevIdx = 0;
				uint32_t nextIdx = 0;

				if(m_jump.find(row) == m_jump.end()) {
					prevIdx = row;
					curIdx = m_jump[row] = row;
				} else {
					prevIdx = m_jump[row];
					curIdx = m_jump[row] = m_nextJump++;
				}
				
				auto it0 = m_cache[row].begin();
				auto it1 = it0;
				std::advance(it1, count);
				std::list<LASPoint*> pts(it0, it1);
				m_cache[row].erase(it0, it1);

				if(prevIdx != curIdx) {
					std::fseek(m_file, prevIdx *  rowSize + 2 * sizeof(uint32_t), SEEK_SET);
					std::fwrite((const void *) &curIdx, sizeof(uint32_t), 1, m_file);
				}

				std::fseek(m_file, curIdx * rowSize, SEEK_SET);
				std::fwrite((const void *) &row, sizeof(uint32_t), 1, m_file);
				std::fwrite((const void *) &count, sizeof(uint32_t), 1, m_file);
				std::fwrite((const void *) &nextIdx, sizeof(uint32_t), 1, m_file);

				for(const LASPoint *pt : pts) {
					pt->write(m_file);
					delete pt;
				}
			}
		}
	}
}

void SortedPointStream::init() {

	if(m_inited) return;
	m_inited = true;

	m_bounds.collapse();

	if(m_rebuild) {

		liblas::ReaderFactory rf;
		std::queue<std::string> m_fileq;
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
				m_fileq.push(file);
			}
		}
		m_bounds.snap(g_abs(m_blockSize));

		m_rowCount = (uint32_t) (m_bounds.height() / m_blockSize) + 1;
		m_running = true;
		m_finalize = false;
		m_filename = Util::tmpFile("/tmp");
		m_nextJump = m_rowCount;

		m_file = std::fopen(m_filename.c_str(), "wb");
		if(!m_file)
			g_runerr("Failed to open point file.");

		std::thread consumer(_consume, this);

		std::list<std::thread> producers;
		for(uint32_t i = 0; i < 1; ++i)
			producers.push_back(std::thread(_produce, this));

		consumer.join();

		m_finalize = true;

		for(std::thread &t : producers)
			t.join();

		std::fclose(m_file);

	}
}

const Bounds& SortedPointStream::bounds() const {
	return m_bounds;
}

bool SortedPointStream::next(std::list<std::shared_ptr<LASPoint> > &pts) {
	if(!m_file) {
		m_file = std::fopen(m_filename.c_str(), "rb");
		m_row = 0;
		if(!m_file) 
			g_argerr("Failed to open " << m_filename);
	}

	uint32_t rowSize = 3 * sizeof(uint32_t) + 1024 * LASPoint::dataSize;
	uint32_t row, count, jump;

	jump = m_row;
	do {
		std::fseek(m_file, jump * rowSize, SEEK_SET);
		std::fread((void *) &row, sizeof(uint32_t), 1, m_file);
		std::fread((void *) &count, sizeof(uint32_t), 1, m_file);
		std::fread((void *) &jump, sizeof(uint32_t), 1, m_file);
		for(uint32_t i = 0; i < count; ++i) {
			std::shared_ptr<LASPoint> pt(new LASPoint());
			pt->read(m_file);
			pts.push_back(pt);
		}
	} while(jump > 0);

	return pts.size() > 0;
}
