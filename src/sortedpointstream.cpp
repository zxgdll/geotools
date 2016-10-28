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

#define CACHE_LEN 4096

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
	void *buf = calloc(LASPoint::dataSize, 1);
	std::fread(buf, LASPoint::dataSize, 1, str);
	read(buf);
	std::free(buf);
}

void LASPoint::write(std::FILE *str) {
	void *buf = calloc(LASPoint::dataSize, 1);
	write(buf);
	std::fwrite(buf, LASPoint::dataSize, 1, str);
	std::free(buf);
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

LASPoint::~LASPoint() {
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

uint32_t _row(const LASPoint &pt, const Bounds &bounds, double blockSize) {
	if(blockSize < 0) {
		return (uint32_t) ((pt.y - bounds.maxy()) / blockSize);
	} else {
		return (uint32_t) ((pt.y - bounds.miny()) / blockSize);
	}
}

void SortedPointStream::produce() {

	liblas::ReaderFactory rf;
	while(true) {

		// Get the file. If the queue is empty, quit.
		std::string file;
		m_fmtx.lock();
		if(m_fileq.size()) {
			file = m_fileq.front();
			m_fileq.pop();
		}
		m_fmtx.unlock();

		if(file.empty()) 
			break;

		// Prepare the LAS source.
		g_debug(" -- reading file " << file);
		std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
		liblas::Reader lasReader = rf.CreateWithStream(instr);
		liblas::Header lasHeader = lasReader.GetHeader();

		// Read each point and add to the cache for its row.
		uint32_t row;
		LASPoint pt;
		while(lasReader.ReadNextPoint()) {
			pt.update(lasReader.GetPoint());
			row = _row(pt, m_bounds, m_blockSize);
			m_cmtx.lock();
			m_cache[row].push_back(new LASPoint(pt));
			m_cmtx.unlock();
		}
	}

}

void SortedPointStream::consume() {

	uint64_t bufLen = 3 * sizeof(uint32_t) + CACHE_LEN * LASPoint::dataSize;
	void *buf = calloc(bufLen, 1);
	std::list<LASPoint*> pts;

	while(m_running || m_cache.size()) {

		uint32_t row;
		uint32_t count = 0;

		m_cmtx.lock();
		for(const auto &it : m_cache) {
			if(!m_running || it.second.size() >= CACHE_LEN) {
				row   = it.first;
				count = g_min(CACHE_LEN, it.second.size());
				auto it0 = m_cache[row].begin();
				pts.resize(0);
				pts.splice(pts.begin(), m_cache[row], it0, std::next(it0, count));
				break;
			}
		}
		m_cmtx.unlock();

		if(!count) {
			m_cmtx.lock();
			m_cache.erase(row);
			m_cmtx.unlock();
			continue;
		}

		uint32_t curIdx  = 0;
		uint32_t prevIdx = 0;
		uint32_t nextIdx = 0;

		m_cmtx.lock();
		if(m_jump.find(row) == m_jump.end()) {
			prevIdx = row;
			curIdx = m_jump[row] = row;
		} else {
			prevIdx = m_jump[row];
			curIdx = m_jump[row] = m_nextJump++;
		}
		m_cmtx.unlock();
	
		char *buf0 = (char *) buf;

		*((uint32_t *) buf0) = row;     buf0 += sizeof(uint32_t);
		*((uint32_t *) buf0) = count;   buf0 += sizeof(uint32_t);
		*((uint32_t *) buf0) = nextIdx; buf0 += sizeof(uint32_t);

		g_debug(" -- writing row " << row << "; " << pts.size() << "; " << m_cache[row].size() << "; " << std::this_thread::get_id());
		for(const LASPoint *pt : pts) {
			pt->write((void *) buf0); buf0 += LASPoint::dataSize;
			delete pt;
		}

		m_wmtx.lock();
		if(prevIdx != curIdx) {
			std::fseek(m_file, prevIdx *  bufLen + 2 * sizeof(uint32_t), SEEK_SET);
			std::fwrite((const void *) &curIdx, sizeof(uint32_t), 1, m_file);
		}
		std::fseek(m_file, curIdx * bufLen, SEEK_SET);
		std::fwrite(buf, bufLen, 1, m_file);
		m_wmtx.unlock();

	}

	free(buf);
}

void SortedPointStream::init() {

	if(m_inited) return;
	m_inited = true;

	m_bounds.collapse();

	if(m_rebuild) {

		liblas::ReaderFactory rf;
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

		uint32_t rowSize = 3 * sizeof(uint32_t) + CACHE_LEN * LASPoint::dataSize;

		m_rowCount = (uint32_t) (m_bounds.height() / g_abs(m_blockSize)) + 1;
		m_running = true;
		m_filename = Util::tmpFile("/tmp");
		m_nextJump = m_rowCount;

		m_file = std::fopen(m_filename.c_str(), "wb");
		if(!m_file)
			g_runerr("Failed to open point file.");
		void *buf = std::calloc(rowSize, 1);
		for(uint32_t i = 0; i < m_rowCount; ++i)
			std::fwrite((const void *) buf, rowSize, 1, m_file);
		free(buf);

		std::list<std::thread> consumers;
		for(uint32_t i = 0; i < 3; ++i)
			consumers.push_back(std::thread(&SortedPointStream::consume, this));

		std::list<std::thread> producers;
		for(uint32_t i = 0; i < 3; ++i)
			producers.push_back(std::thread(&SortedPointStream::produce, this));

		for(std::thread &t : producers)
			t.join();

		m_running = false;

		for(std::thread &t : consumers)
			t.join();

		g_debug(" -- sorting done.");

		std::fclose(m_file);
		m_file = nullptr;

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

	uint32_t rowSize = 3 * sizeof(uint32_t) + CACHE_LEN * LASPoint::dataSize;
	uint32_t row, count, jump;

	jump = m_row;
	do {
		if(std::fseek(m_file, jump * rowSize, SEEK_SET))
			return false;
		// TODO: Read to a buffer		
		std::fread((void *) &row, sizeof(uint32_t), 1, m_file);
		std::fread((void *) &count, sizeof(uint32_t), 1, m_file);
		std::fread((void *) &jump, sizeof(uint32_t), 1, m_file);
		for(uint32_t i = 0; i < count; ++i) {
			std::shared_ptr<LASPoint> pt(new LASPoint());
			pt->read(m_file);
			pts.push_back(pt);
		}
	} while(jump > 0);

	++m_row;

	return pts.size() > 0;
}
