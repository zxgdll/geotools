#ifndef __LASREADER_HPP__
#define __LASREADER_HPP__

#include <cstdio>
#include <vector>
#include <memory>

#include "util.hpp"
#include "laspoint.hpp"

using namespace geotools::las;
using namespace geotools::util;

class LASReader {
private:
	std::FILE *m_f;
	std::string m_file;
	uint16_t m_sourceId;
	std::string m_version;
	uint16_t m_headerSize;
	uint32_t m_offset;
	uint8_t m_pointFormat;
	uint16_t m_pointLength;
	uint64_t m_pointCount;
	uint64_t m_pointCountByReturn[5];
	double m_xScale;
	double m_yScale;
	double m_zScale;
	double m_xOffset;
	double m_yOffset;
	double m_zOffset;
	double m_xMin;
	double m_yMin;
	double m_zMin;
	double m_xMax;
	double m_yMax;
	double m_zMax;

	uint64_t m_curPoint;

	std::unique_ptr<Buffer> m_buf;
	std::queue<LASPoint> m_pts;

	void load() {
		m_f = std::fopen(m_file.c_str(), "rb");
		if(!m_f)
			g_runerr("Failed to open " << m_file);

		std::fseek(m_f, 4, SEEK_SET);
		std::fread((void *) &m_sourceId, sizeof(uint16_t), 1, m_f);

		std::fseek(m_f, 24, SEEK_SET);
		char maj, min;
		std::fread((void *) &maj, sizeof(char), 1, m_f);
		std::fread((void *) &min, sizeof(char), 1, m_f);
		//m_version = "" + maj + "." + min;

		std::fseek(m_f, 94, SEEK_SET);
		std::fread((void *) &m_headerSize, sizeof(uint16_t), 1, m_f);
	
		std::fseek(m_f, 96, SEEK_SET);
		std::fread((void *) &m_offset, sizeof(uint32_t), 1, m_f);

		std::fseek(m_f, 104, SEEK_SET);
		std::fread((void *) &m_pointFormat, sizeof(uint8_t), 1, m_f);
		std::fread((void *) &m_pointLength, sizeof(uint16_t), 1, m_f);

		uint32_t pc, pcbr[5];
		std::fread((void *) &pc, sizeof(uint32_t), 1, m_f);
		std::fread((void *) &pcbr, 5 * sizeof(uint32_t), 1, m_f);
		m_pointCount = pc;
		for(int i = 0; i < 5; ++i)
			m_pointCountByReturn[i] = pcbr[i];

		std::fread((void *) &m_xScale, sizeof(double), 1, m_f);
		std::fread((void *) &m_yScale, sizeof(double), 1, m_f);
		std::fread((void *) &m_zScale, sizeof(double), 1, m_f);

		std::fread((void *) &m_xOffset, sizeof(double), 1, m_f);
		std::fread((void *) &m_yOffset, sizeof(double), 1, m_f);
		std::fread((void *) &m_zOffset, sizeof(double), 1, m_f);

		std::fread((void *) &m_xMax, sizeof(double), 1, m_f);
		std::fread((void *) &m_xMin, sizeof(double), 1, m_f);
		std::fread((void *) &m_yMax, sizeof(double), 1, m_f);
		std::fread((void *) &m_yMin, sizeof(double), 1, m_f);
		std::fread((void *) &m_zMax, sizeof(double), 1, m_f);
		std::fread((void *) &m_zMin, sizeof(double), 1, m_f);

		// TODO: Extended point count.

		LASPoint::setScale(m_xScale, m_yScale, m_zScale);
		reset();
	}

public:
	LASReader(const std::string &file) :
		m_f(nullptr),
		m_file(file) {
		load();
	}

	~LASReader() {
		if(m_f)
			std::fclose(m_f);
		m_buf.release();
	}

	size_t m_batchSize;
	size_t BATCH_SIZE = 100000;

	void reset() {
		m_curPoint = 0;
		std::fseek(m_f, m_offset, SEEK_SET);
	}

	bool loadBatch() {
		if(!m_f || m_curPoint >= m_pointCount) 
			return false;
		m_batchSize = g_min(BATCH_SIZE, m_pointCount - m_curPoint);
		m_buf.reset(new Buffer(m_batchSize * m_pointLength));
		g_debug(" -- loading " << m_batchSize << " points");
		std::fread((void *) m_buf->buf, m_batchSize * m_pointLength, 1, m_f);
		return true;
	}

	bool next(LASPoint &pt) {
		if(m_curPoint >= m_pointCount)
			return false;
		if(m_curPoint % BATCH_SIZE == 0 && !loadBatch())
			return false;
		char *buf = ((char *) m_buf->buf) + (m_curPoint % BATCH_SIZE) * m_pointLength;
		pt.readLAS(buf, m_pointFormat);
		++m_curPoint;
		return true;
	}

	Bounds bounds() {
		return Bounds(m_xMin, m_yMin, m_xMax, m_yMax, m_zMin, m_zMax);
	}

	uint64_t pointCount() {
		return m_pointCount;
	}

};

#endif
