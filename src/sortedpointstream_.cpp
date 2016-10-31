#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <memory>
#include <cstdio>
#include <cerrno>

#include "liblas/liblas.hpp"

#include "geotools.h"
#include "sortedpointstream.hpp"
#include "lasutil.hpp"

// The number of points to store in a row before writing.
// TODO: Configure/optimize.
#define CACHE_LEN 2048

// The size in bytes of the header row.
#define HEADER_LEN 64
 
using namespace geotools::las;
using namespace geotools::util;

double LASPoint::scaleX = 0;
double LASPoint::scaleZ = 0;
double LASPoint::scaleY = 0;

// Return the row number, given the point's y position,
// bounds and block size.
uint64_t _row(const std::shared_ptr<LASPoint> &pt, const Bounds &bounds, double blockSize) {
	if(blockSize < 0) {
		return (uint64_t) ((pt->y - bounds.maxy()) / blockSize);
	} else {
		return (uint64_t) ((pt->y - bounds.miny()) / blockSize);
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
	returnNum = pt.GetReturnNumber();
	numReturns = pt.GetNumberOfReturns();
	cls = pt.GetClassification().GetClass();
	scanAngle = pt.GetScanAngleRank();
}

void LASPoint::read(std::istream &str) {
	int32_t xx, yy, zz;
	str.read((char *) &xx, sizeof(int32_t));
	str.read((char *) &yy, sizeof(int32_t));
	str.read((char *) &zz, sizeof(int32_t));
	str.read((char *) &intensity, sizeof(uint16_t));
	str.read((char *) &returnNum, sizeof(uint16_t));
	str.read((char *) &numReturns, sizeof(uint16_t));
	str.read((char *) &cls, sizeof(uint8_t));
	str.read((char *) &scanAngle, sizeof(int8_t));
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
	str.write((char *) &returnNum, sizeof(uint16_t));
	str.write((char *) &numReturns, sizeof(uint16_t));
	str.write((char *) &cls, sizeof(uint8_t));
	str.write((char *) &scanAngle, sizeof(int8_t));
}

void LASPoint::read(std::FILE *str) {
	Buffer buffer(LASPoint::dataSize);
	std::fread(buffer.buf, LASPoint::dataSize, 1, str);
	read(buffer.buf);
}

void LASPoint::write(std::FILE *str) {
	Buffer buffer(LASPoint::dataSize);
	write(buffer.buf);
	std::fwrite(buffer.buf, LASPoint::dataSize, 1, str);
}

void LASPoint::read(void *str) {
	char *ptr  = (char *) str;
	x          = *((int32_t *)  ptr) * scaleX; ptr += sizeof(int32_t);
	y          = *((int32_t *)  ptr) * scaleY; ptr += sizeof(int32_t);
	z          = *((int32_t *)  ptr) * scaleZ; ptr += sizeof(int32_t);
	intensity  = *((uint16_t *) ptr);          ptr += sizeof(uint16_t);
	returnNum  = *((uint16_t *) ptr);          ptr += sizeof(uint16_t);
	numReturns = *((uint16_t *) ptr);          ptr += sizeof(uint16_t);
	cls        = *((uint8_t *)  ptr);          ptr += sizeof(uint8_t);
	scanAngle  = *((int8_t *)   ptr);
}

void LASPoint::write(void *str) const {
	char *ptr = (char *) str;
	*((int32_t *)  ptr) = (int32_t) (x / scaleX); ptr += sizeof(int32_t);
	*((int32_t *)  ptr) = (int32_t) (y / scaleY); ptr += sizeof(int32_t);
	*((int32_t *)  ptr) = (int32_t) (z / scaleZ); ptr += sizeof(int32_t);
	*((uint16_t *) ptr) = intensity;              ptr += sizeof(uint16_t);
	*((uint16_t *) ptr) = returnNum;              ptr += sizeof(uint16_t);
	*((uint16_t *) ptr) = numReturns;             ptr += sizeof(uint16_t);
	*((uint8_t *)  ptr) = cls;                    ptr += sizeof(uint8_t);
	*((int8_t *)   ptr) = scanAngle;
}

bool LASPoint::last() const {
	return numReturns > 0 && returnNum == numReturns;
}

bool LASPoint::first() const {
	return numReturns > 0 && returnNum == 1;
}

bool LASPoint::intermediate() const {
	return numReturns > 2 && returnNum > 1 && returnNum < numReturns;
}

bool LASPoint::ground() const {
	return cls == 2;
}

bool LASPoint::single() const {
	return numReturns == 1;
}

LASPoint::~LASPoint() {
}

SortedPointStream::SortedPointStream(const std::list<std::string> &files, const std::string &cacheFile, 
	double blockSize, bool rebuild, bool snap, uint32_t threads) :
	m_files(files),
	m_cacheFile(cacheFile),
	m_blockSize(blockSize),
	m_rebuild(rebuild),
	m_inited(false),
	m_file(nullptr),
	m_snap(snap), 
	m_threads(threads) {
}

SortedPointStream::~SortedPointStream() {
	if(m_file)
		std::fclose(m_file);
}

uint64_t SortedPointStream::pointCount() const {
	return m_pointCount;
}

uint64_t SortedPointStream::rowCount() const {
	return m_rowCount;
}


void SortedPointStream::produce() {

	liblas::ReaderFactory rf;
	while(true) {

		// Get the file. If the queue is empty, quit.
		std::string file;
		m_fmtx.lock();
		if(!m_fileq.size()) {
			m_fmtx.unlock();
			break;
		}
		file = m_fileq.front();
		m_fileq.pop();
		m_fmtx.unlock();

		// Prepare the LAS source.
		g_debug(" -- reading file " << file);
		std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
		liblas::Reader lasReader = rf.CreateWithStream(instr);
		liblas::Header lasHeader = lasReader.GetHeader();

		// Read each point and add to the cache for its row.
		uint64_t row;
		while(lasReader.ReadNextPoint()) {
			std::shared_ptr<LASPoint> pt(new LASPoint());
			pt->update(lasReader.GetPoint());
			row = _row(pt, m_bounds, m_blockSize);
			//g_debug(" -- point for row " << row);
			m_cmtx.lock();
			m_cache[row].push_back(std::move(pt));
			m_flush.push(row);
			m_cmtx.unlock();
		}
	}
}

void SortedPointStream::consume() {

	uint64_t bufLen = 3 * sizeof(uint64_t) + CACHE_LEN * LASPoint::dataSize;
	Buffer buffer(bufLen);
	std::list<std::shared_ptr<LASPoint> > pts;

	while(m_running || m_flush.size()) {

		uint64_t row;

		m_cmtx.lock();
		if(m_flush.size()) {
			row = m_flush.front();
			m_flush.pop();
		} else {
			m_cmtx.unlock();
			continue;
		}
		m_cmtx.unlock();

		uint64_t count = 0;
		uint64_t curIdx  = 0;
		uint64_t prevIdx = 0;
		uint64_t nextIdx = 0;

		m_cmtx.lock();
		if(m_cache[row].size() >= CACHE_LEN || !m_running) {

			count = g_min(CACHE_LEN, m_cache[row].size());
			auto it0 = m_cache[row].begin();
			pts.resize(0);
			pts.splice(pts.begin(), m_cache[row], it0, std::next(it0, count));
			
			if(m_cache[row].size() >= CACHE_LEN || (!m_running && m_cache[row].size() > 0)) {
				//g_debug(" -- " << m_cache[row].size() << " " << m_running);
				for(uint32_t i = 0; i < (m_cache[row].size() / CACHE_LEN) + 1; ++i)
					m_flush.push(row);
			}

			if(m_jump.find(row) == m_jump.end()) {
				prevIdx = row;
				curIdx = m_jump[row] = row;
			} else {
				prevIdx = m_jump[row];
				curIdx = m_jump[row] = m_nextJump++;
			}
		}
		m_cmtx.unlock();

		if(!count)
			continue;

		char *buf0 = (char *) buffer.buf;

		*((uint64_t *) buf0) = row;     buf0 += sizeof(uint64_t);
		*((uint64_t *) buf0) = count;   buf0 += sizeof(uint64_t);
		*((uint64_t *) buf0) = nextIdx; buf0 += sizeof(uint64_t);

		//g_debug(" -- writing row " << row << "; " << pts.size() << "; " << m_cache[row].size() << "; " << std::this_thread::get_id());
		for(const std::shared_ptr<LASPoint> &pt : pts) {
			pt->write((void *) buf0); 
			buf0 += LASPoint::dataSize;
		}

		m_wmtx.lock();
		if(prevIdx != curIdx) {
			std::fseek(m_file, prevIdx *  bufLen + HEADER_LEN + 2 * sizeof(uint64_t), SEEK_SET);
			std::fwrite((const void *) &curIdx, sizeof(uint64_t), 1, m_file);
		}
		std::fseek(m_file, curIdx * bufLen + HEADER_LEN, SEEK_SET);
		std::fwrite(buffer.buf, bufLen, 1, m_file);
		m_wmtx.unlock();
	}
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
		if(m_snap)
			m_bounds.snap(g_abs(m_blockSize));

		uint64_t bufSize = 3 * sizeof(uint64_t) + CACHE_LEN * LASPoint::dataSize;

		m_rowCount = (uint64_t) (m_bounds.height() / g_abs(m_blockSize)) + 1;
		m_running = true;
		m_nextJump = m_rowCount;

		// Open file and write all zeroes for all rows.
		m_file = std::fopen(m_cacheFile.c_str(), "wb");
		if(!m_file)
			g_runerr("Failed to open point file.");
		{
			Buffer buffer(bufSize);
			std::fwrite((const void *) buffer.buf, HEADER_LEN, 1, m_file);
			for(uint64_t i = 0; i < m_rowCount; ++i)
				std::fwrite((const void *) buffer.buf, bufSize, 1, m_file);
		}
		
		// Start up a set of consumers to sort points.
		std::list<std::thread> consumers;
		for(uint32_t i = 0; i < 2; ++i)
			consumers.push_back(std::thread(&SortedPointStream::consume, this));

		// Start up a set of consumers to read las files.
		std::list<std::thread> producers;
		for(uint32_t i = 0; i < 1; ++i)
			producers.push_back(std::thread(&SortedPointStream::produce, this));

		// Wait for producers to finish
		for(std::thread &t : producers)
			t.join();

		m_cmtx.lock();
		for(const auto it : m_cache)
			m_flush.push(it.first);
		m_cmtx.unlock();

		// Tell the consumers to finish and wait.
		m_running = false;
		for(std::thread &t : consumers)
			t.join();

		double minx = m_bounds.minx(), miny = m_bounds.miny();
		double maxx = m_bounds.maxx(), maxy = m_bounds.maxy();
		double scaleX = LASPoint::scaleX, scaleY = LASPoint::scaleY, scaleZ = LASPoint::scaleZ;
		// Write the bounds and row count toe the first line in the file.
		std::fseek(m_file, 0, SEEK_SET);
		std::fwrite((const void *) &minx, sizeof(double), 1, m_file);
		std::fwrite((const void *) &miny, sizeof(double), 1, m_file);
		std::fwrite((const void *) &maxx, sizeof(double), 1, m_file);
		std::fwrite((const void *) &maxy, sizeof(double), 1, m_file);
		std::fwrite((const void *) &scaleX, sizeof(double), 1, m_file);
		std::fwrite((const void *) &scaleY, sizeof(double), 1, m_file);
		std::fwrite((const void *) &scaleZ, sizeof(double), 1, m_file);
		std::fwrite((const void *) &m_rowCount, sizeof(uint64_t), 1, m_file);
		
		std::fclose(m_file);
		m_file = nullptr;
		g_debug(" -- sorting done.");
	} else {
		m_file = std::fopen(m_cacheFile.c_str(), "rb");
		if(!m_file)
			g_runerr("Couldn't open cache file for reading.");
		double minx, miny, maxx, maxy, scaleX, scaleY, scaleZ;
		std::fread((void *) &minx, sizeof(double), 1, m_file);
		std::fread((void *) &miny, sizeof(double), 1, m_file);
		std::fread((void *) &maxx, sizeof(double), 1, m_file);
		std::fread((void *) &maxy, sizeof(double), 1, m_file);
		std::fread((void *) &scaleX, sizeof(double), 1, m_file);
		std::fread((void *) &scaleY, sizeof(double), 1, m_file);
		std::fread((void *) &scaleZ, sizeof(double), 1, m_file);
		std::fread((void *) &m_rowCount, sizeof(uint64_t), 1, m_file);
		std::fclose(m_file);
		m_file = nullptr;
		m_bounds.extend(minx, miny);
		m_bounds.extend(maxx, maxy);
		LASPoint::setScale(scaleX, scaleY, scaleZ);
	}
}

const Bounds& SortedPointStream::bounds() const {
	return m_bounds;
}

uint64_t SortedPointStream::bufferSize() const {
	return 3 * sizeof(uint64_t) + CACHE_LEN * LASPoint::dataSize;
}

bool SortedPointStream::next(std::list<std::shared_ptr<LASPoint> > &pts) {
	if(m_mmap.get() == nullptr) {
		g_debug(" -- next opening file");
		m_mmap.reset(Util::mapFile(m_cacheFile, 0).release());
		m_row = 0;
	}

	uint64_t row, count, jump;
	uint64_t bufSize = bufferSize();

	jump = m_row;
	do {
		char *buf = ((char *) m_mmap->data()) + jump * bufSize + HEADER_LEN;
		row   = *((uint64_t *) buf);  buf += sizeof(uint64_t);
		count = *((uint64_t *) buf);  buf += sizeof(uint64_t);
		jump  = *((uint64_t *) buf);  buf += sizeof(uint64_t);
		if(row != m_row && (count > 0 || jump > 0))
			g_runerr("Invalid row header " << row << "; " << m_row << " expected.");

		for(uint64_t i = 0; i < count; ++i) {
			std::shared_ptr<LASPoint> pt(new LASPoint());
			pt->read(buf);
			pts.push_back(std::move(pt));
			buf += LASPoint::dataSize;
		}
	} while(jump > 0);

	return ++m_row < m_rowCount;
}
