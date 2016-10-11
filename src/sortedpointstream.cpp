#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <memory>
#include <cstdio>

#include "geos/geom/Geometry.h"
#include "geos/geom/Polygon.h"
#include "geos/geom/Coordinate.h"

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

SortedPointStream::SortedPointStream(const std::string &file, unsigned int numBlocks) :
	m_numBlocks(numBlocks), 
	m_file(file),
	m_currentBlock(0),
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

unsigned int SortedPointStream::pointCount() {
	if(!m_inited)
		g_runerr("Not inited.");
	return m_pointCount;
}

void SortedPointStream::init() {
	#pragma omp critical(__sps_init)
	{
		if(!m_inited) {
			g_debug(" -- init: opening file: " << m_file);
			std::ifstream instr(m_file.c_str(), std::ios::in|std::ios::binary);
			liblas::ReaderFactory rf;
			liblas::Reader lasReader = rf.CreateWithStream(instr);
			liblas::Header lasHeader = lasReader.GetHeader();

			LASPoint::setScale(lasHeader.GetScaleX(), lasHeader.GetScaleY(), lasHeader.GetScaleZ());
			m_pointCount = lasHeader.GetPointRecordsCount();

			g_debug(" -- init computing bounds");
			if(!LasUtil::computeLasBounds(lasHeader, m_fileBounds, 2))
				LasUtil::computeLasBounds(lasReader, m_fileBounds, 2); // If the header bounds are bogus.

			std::map<unsigned int, std::unique_ptr<std::ofstream> > blockFiles;
			std::map<unsigned int, std::string> blockFilenames;

			g_debug(" -- init opening out files");
			for(unsigned int i = 0; i <= m_numBlocks; ++i) { // N.B. Extra block at the end.
				std::stringstream fs;
				fs << "/tmp/sps_" << i;
				std::string file = fs.str();
				std::unique_ptr<std::ofstream> os(new std::ofstream(file, std::ios::out|std::ios::binary));
				blockFiles[i] = std::move(os);
				m_blockFilenames[i] = file;
			}

			LASPoint pt;
			unsigned int block;

			g_debug(" -- init writing points");
			int cnt = 0;
			while(lasReader.ReadNextPoint()) {
				++cnt;
				pt.update(lasReader.GetPoint());
				block = (unsigned int) ((pt.y - m_fileBounds.miny()) / m_fileBounds.height() * m_numBlocks);
				pt.write(*(blockFiles[block].get()));
				m_cacheCounts[block]++;
			}

			for(const auto &it : m_cacheCounts)
				g_debug(" -- init block " << it.first << "; count " << it.second);

			g_debug(" -- init opening in files");
			for(unsigned int i = 0; i < m_numBlocks; ++i) {
				blockFiles[i]->close();
				std::unique_ptr<std::ifstream> is(new std::ifstream(m_blockFilenames[i], std::ios::in|std::ios::binary));
				m_cacheFiles[i] = std::move(is);
			}

			m_currentBlock = 0;
			m_inited = true;
		}
	} // omp
}

Bounds SortedPointStream::fileBounds() {
	if(!m_inited)
		g_runerr("Not inited.");
	return m_fileBounds;
}

Bounds SortedPointStream::currentBounds() {
	if(!m_inited)
		g_runerr("Not inited.");
	return m_fileBounds;
}

bool SortedPointStream::next(LASPoint &pt, geos::geom::Geometry **geom) {

	if(!m_inited)
		g_runerr("Not inited.");

	if(m_currentBlock >= m_numBlocks) 
		return false;

	if(m_cacheCounts[m_currentBlock] == 0) {
		g_debug(" -- next - switching to next block " << m_currentBlock);
		geos::geom::Geometry *poly = (geos::geom::Geometry *) SimpleGeom::createPolygon(m_fileBounds);
		if(*geom == nullptr) {
			*geom = poly;
		} else {
			geos::geom::Geometry *tmp;
			if((*geom)->intersects(poly)) {
				tmp = *geom;
				*geom = (*geom)->intersection(poly);
				delete tmp;
			} else {
				tmp = *geom;
				*geom = (*geom)->Union(poly);
				delete geom;
			}
		}
		m_currentBlock++;
	}

	if(m_currentBlock >= m_numBlocks) 
		return false;

	pt.read(*(m_cacheFiles[m_currentBlock].get()));
	m_cacheCounts[m_currentBlock]--;

	return true;
}
