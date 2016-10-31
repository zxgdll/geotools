#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <memory>
#include <cstdio>
#include <cerrno>
#include <parallel/algorithm>

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

	liblas::ReaderFactory rf;
	Bounds workBounds; // Configure
	for(const std::string &file : m_files) {
		g_debug(" -- init: opening file: " << file);
		std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
		liblas::Reader lasReader = rf.CreateWithStream(instr);
		liblas::Header lasHeader = lasReader.GetHeader();

		g_debug(" -- init computing bounds");
		Bounds fileBounds;
		if(!LasUtil::computeLasBounds(lasHeader, fileBounds, 2))
			LasUtil::computeLasBounds(lasReader, fileBounds, 2); // If the header bounds are bogus.
		if(fileBounds.intersects(workBounds)) {
			m_bounds.extend(fileBounds);
			m_pointCount += lasHeader.GetPointRecordsCount();
		}
	}
	if(m_snap)
		m_bounds.snap(g_abs(m_blockSize));

	m_rowCount = (uint64_t) (m_bounds.height() / g_abs(m_blockSize)) + 1;
	m_row = 0;
	m_idx = 0;

	size_t s = sizeof(std::vector<LASPoint>) + sizeof(LASPoint) * m_pointCount * 1.1;
	m_mfile.reset(Util::mapFile("cache.tmp", s).release());
	void *ptr = m_mfile->data();
	m_lst = (LASPoint*) ptr;

	try {
		size_t idx = 0;
		std::vector<std::string> files(m_files.begin(), m_files.end());
		#pragma omp parallel
		{
		
			#pragma omp for
			for(uint32_t i = 0; i < files.size(); ++i) {
				const std::string &file = files[i];
				g_debug(" -- init: opening file: " << file);
				std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
				liblas::Reader lasReader = rf.CreateWithStream(instr);
				liblas::Header lasHeader = lasReader.GetHeader();
				while(lasReader.ReadNextPoint()) {
					LASPoint pt;
					pt.update(lasReader.GetPoint());
					m_lst[idx] = std::move(pt);
					#pragma omp atomic
					idx++;
				}
			}
		}
	} catch(const std::bad_alloc &c) {
		g_runerr(" -- bad_alloc " << c.what());
	}

	struct {
		bool operator()(const LASPoint &p1, const LASPoint &p2) {
			return p1.y > p2.y; // TODO: THis is because of negative res-y
		}
	} sortfn;
	g_debug("sorting");
	std::sort(m_lst, m_lst + m_pointCount, sortfn);
}

const Bounds& SortedPointStream::bounds() const {
	return m_bounds;
}

uint64_t SortedPointStream::bufferSize() const {
	return 3 * sizeof(uint64_t) + CACHE_LEN * LASPoint::dataSize;
}

bool SortedPointStream::next(std::list<std::shared_ptr<LASPoint> > &pts) {
	if(m_row >= m_rowCount)
		return false;
	for(uint64_t i = m_idx; i < m_pointCount; ++i) {
		std::shared_ptr<LASPoint> pt(new LASPoint(*(m_lst + i)));
		uint32_t row = _row(pt, m_bounds, m_blockSize);
		if(row == m_row) {
			pts.push_back(std::move(pt));
		} else {
			++m_row;
			m_idx = i;
			return true;
		}
	}
	++m_row; // This happens at the end of the list.
	return false;
}
