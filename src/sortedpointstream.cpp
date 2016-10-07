#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <memory>

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
	// close streams
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
				blockFilenames[i] = file;
			}

			LASPoint pt;
			unsigned int block;

			g_debug(" -- init writing points");
			while(lasReader.ReadNextPoint()) {
				pt.update(lasReader.GetPoint());
				block = (unsigned int) ((pt.y - m_fileBounds.miny()) / m_fileBounds.height() * m_numBlocks);
				pt.write(*(blockFiles[block].get()));
				m_cacheCounts[block]++;
			}

			g_debug(" -- init opening in files");
			for(unsigned int i = 0; i < m_numBlocks; ++i) {
				blockFiles[i]->close();
				std::unique_ptr<std::ifstream> is(new std::ifstream(blockFilenames[i], std::ios::in|std::ios::binary));
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
		std::unique_ptr<geos::geom::Polygon> poly = std::move(SimpleGeom::createPolygon(m_fileBounds));
		if(*geom == nullptr) {
			*geom = (geos::geom::Geometry *) poly.release();
		} else {
			const geos::geom::Geometry *g = (geos::geom::Geometry *) poly.release();
			if((*geom)->intersects(g)) {
				*geom = (*geom)->intersection(g);
			} else {
				*geom = (*geom)->Union(g);
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
