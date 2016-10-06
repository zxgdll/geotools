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

SortedPointStream::SortedPointStream(const std::string &file, unsigned int numBlocks, 
		unsigned char mode, bool mem) :
	m_numBlocks(numBlocks), 
	m_mode(mode), 
	m_mem(mem),
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

			m_pointCount = lasHeader.GetPointRecordsCount();

			g_debug(" -- init computing bounds");
			if(!LasUtil::computeLasBounds(lasHeader, m_fileBounds, 2))
				LasUtil::computeLasBounds(lasReader, m_fileBounds, 2); // If the header bounds are bogus.

			double blockHeight = m_fileBounds.height() / m_numBlocks;
			std::map<unsigned int, std::unique_ptr<std::ofstream> > blockFiles;
			std::map<unsigned int, std::string> blockFilenames;

			g_debug(" -- init opening out files");
			for(unsigned int i = 0; i < m_numBlocks; ++i) {
				std::stringstream fs;
				fs << "/tmp/sps_" << i;
				std::string file = fs.str();
				std::unique_ptr<std::ofstream> os(new std::ofstream(file, std::ios::out|std::ios::binary));
				blockFiles[i] = std::move(os);
				blockFilenames[i] = file;
			}

			g_debug(" -- init writing points");
			while(lasReader.ReadNextPoint()) {
				liblas::Point pt = lasReader.GetPoint();

				// Get the point coordinates.
				double px = pt.GetX();
				double py = pt.GetY();
				double pz;
				if(m_mode & SortedPointStream::intensity) {
					pz = pt.GetIntensity();
				} else {
					pz = pt.GetZ();
				}

				unsigned int block = (unsigned int) ((py - m_fileBounds.miny()) / m_fileBounds.height() * m_numBlocks);
				*(blockFiles[block].get()) << px << py << pz;
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

geos::geom::GeometryFactory _gf;

bool SortedPointStream::next(Point &pt, geos::geom::Geometry **geom) {
	if(!m_inited)
		g_runerr("Not inited.");
	if(m_currentBlock >= m_numBlocks) 
		return false;
	if(m_cacheCounts[m_currentBlock] == 0)
		m_currentBlock++;
	if(m_currentBlock >= m_numBlocks) 
		return false;

	*(m_cacheFiles[m_currentBlock].get()) >> pt.x >> pt.y >> pt.z;
	m_cacheCounts[m_currentBlock]--;

	if(m_cacheCounts[m_currentBlock] == 0) {
		std::unique_ptr<geos::geom::Polygon> poly = SimpleGeom::createPolygon(m_fileBounds);
		const geos::geom::Geometry *g = (geos::geom::Geometry *) poly.get();
		if((*geom)->intersects(g)) {
			*geom = (*geom)->intersection(g);
		} else {
			*geom = (*geom)->Union(g);
		}
	}
	return true;
}
