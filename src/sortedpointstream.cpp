#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <memory>
#include <cstdio>

#include "liblas/liblas.hpp"

#include "geotools.h"
#include "sortedpointstream.hpp"
#include "simplegeom.hpp"
#include "lasutil.hpp"

using namespace geotools::las;
using namespace geotools::util;
using namespace geotools::geom;

double LASPoint::scaleX = 0;
double LASPoint::scaleZ = 0;
double LASPoint::scaleY = 0;

SortedPointStream::SortedPointStream(const std::list<std::string> &files, double blockSize) :
	m_blockSize(blockSize),
	m_files(files),
	m_inited(false) {
}

SortedPointStream::~SortedPointStream() {
	for(const auto &it : m_blockFilenames) {
		try {
			std::remove(it.second.c_str());
		} catch(...) {
			g_warn("Failed to delete " << it.second);
		}
	}
}

unsigned int SortedPointStream::pointCount() const {
	if(!m_inited)
		g_runerr("Not inited.");
	return m_pointCount;
}

void SortedPointStream::init() {
	#pragma omp critical(__sps_init)
	{
		if(!m_inited) {

			m_bounds.collapse();
			m_pointCount = 0;

			for(const std::string &file : m_files) {
				g_debug(" -- init: opening file: " << file);
				std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
				liblas::ReaderFactory rf;
				liblas::Reader lasReader = rf.CreateWithStream(instr);
				liblas::Header lasHeader = lasReader.GetHeader();

				LASPoint::setScale(lasHeader.GetScaleX(), lasHeader.GetScaleY(), lasHeader.GetScaleZ());

				g_debug(" -- init computing bounds");
				Bounds fileBounds;
				if(!LasUtil::computeLasBounds(lasHeader, fileBounds, 2))
					LasUtil::computeLasBounds(lasReader, fileBounds, 2); // If the header bounds are bogus.
				m_bounds.extend(fileBounds);
			}

			m_minx = m_bounds.minx();
			m_miny = m_bounds.miny();
			m_cols = (int) (m_bounds.width() / m_blockSize) + 1;
			m_rows = (int) (m_bounds.height() / m_blockSize) + 1;
			m_numBlocks = m_cols * m_rows;

			g_debug(" -- init creating out files");
			for(unsigned int row = 0; row <= m_rows; ++row) {
				for(unsigned int col = 0; col <= m_cols; ++col) {
					std::stringstream fs;
					fs << "/tmp/sps_" << row << "_" << col << ".tmp";
					unsigned long block = ((unsigned long) row << 32) | col;
					m_blockFilenames[block] = fs.str();
				}
			}

			LASPoint pt;
			std::map<unsigned long, std::unique_ptr<std::ofstream> > blockFiles;

			for(const std::string &file : m_files) {

				std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
				liblas::ReaderFactory rf;
				liblas::Reader lasReader = rf.CreateWithStream(instr);
				liblas::Header lasHeader = lasReader.GetHeader();

				g_debug(" -- init writing points");
				int cnt = 0;
				while(lasReader.ReadNextPoint()) {
					++cnt;
					pt.update(lasReader.GetPoint());
					unsigned int row = (unsigned int) pt.y / m_rows;
					unsigned int col = (unsigned int) pt.x / m_cols;
					unsigned long block = ((unsigned long) row << 32) | col;
					if(blockFiles.find(block) == blockFiles.end()) {
						std::unique_ptr<std::ofstream> os(new std::ofstream(m_blockFilenames[block], std::ios::out|std::ios::binary));
						blockFiles[block] = std::move(os);
					}
					pt.write(*(blockFiles[block].get()));
					m_cacheCounts[block]++;
					++m_pointCount;
				}
			}

			g_debug(" -- init closing out files");
			for(const auto &it : blockFiles)
				it.second->close();

			m_row = 0;
			m_col = 0;
			m_inited = true;
		}
	} // omp
}

bool SortedPointStream::contains(double x, double y) const {
	return m_boundsTracker.contains(x, y);
}

bool SortedPointStream::contains(double x1, double y1, double x2, double y2) const {
	return m_boundsTracker.contains(x1, y1, x2, y2);
}

const Bounds& SortedPointStream::bounds() const {
	return m_bounds;
}

bool SortedPointStream::next(LASPoint &pt) {

	if(!m_inited)
		g_runerr("Not inited.");

	unsigned long block = ((unsigned long) m_row << 32) | m_col;

	while(m_cacheCounts[block] == 0) {
		g_debug(" -- next - switching to next block " << block);
		m_boundsTracker.add(m_minx + m_col * m_blockSize, m_miny + m_row * m_blockSize, 
			m_minx + (m_col + 1) * m_blockSize, m_miny + (m_row + 1) * m_blockSize);

		++m_col;

		if(m_col >= m_cols) {
			m_col = 0;
			++m_row;
		}

		if(m_row >= m_rows && m_col >= m_cols)
			return false;

	}


	block = ((unsigned long) m_row << 32) | m_col;
	pt.read(*(m_cacheFiles[block].get()));
	m_cacheCounts[block]--;

	return true;
}
