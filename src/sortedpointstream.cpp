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

std::string _rowFile(unsigned int row) {
	std::stringstream fs;
	fs << "/tmp/sps_" << row << ".tmp";
	return fs.str();
}

unsigned int _getRow(double y, const Bounds &bounds, double blockSize) {
	if(blockSize < 0) {
		return (unsigned int) ((y - bounds.maxy()) / blockSize);
	} else {
		return (unsigned int) ((bounds.miny() - y) / blockSize);
	}
}

SortedPointStream::SortedPointStream(const std::list<std::string> &files, double blockSize, bool rebuild) {
	init(files, blockSize, rebuild);
}

SortedPointStream::~SortedPointStream() {
	delete m_pf;
}

unsigned int SortedPointStream::pointCount() const {
	return m_pf->pointCount();
}

unsigned int SortedPointStream::rowCount() const {
	return m_pf->rowCount();
}

void SortedPointStream::init(const std::list<std::string> &files, double blockSize, bool rebuild) {

	m_bounds.collapse();
	m_pointCount = 0;
	m_row = 0;

	liblas::ReaderFactory rf;
	std::unordered_map<unsigned int, std::string> rowFileNames;

	for(const std::string &file : files) {
		g_debug(" -- init: opening file: " << file);
		std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
		liblas::Reader lasReader = rf.CreateWithStream(instr);
		liblas::Header lasHeader = lasReader.GetHeader();

		LASPoint::setScale(lasHeader.GetScaleX(), lasHeader.GetScaleY(), lasHeader.GetScaleZ());

		g_debug(" -- init computing bounds");
		Bounds fileBounds;
		if(!LasUtil::computeLasBounds(lasHeader, fileBounds, 2))
			LasUtil::computeLasBounds(lasReader, fileBounds, 2); // If the header bounds are bogus.
		m_bounds.extend(fileBounds);
	}
	m_bounds.snap(g_abs(blockSize));

	m_pf = new PointFile("/tmp/pf.tmp", m_bounds, blockSize);
	m_pf->openWrite();

	const std::vector<std::string> files0(files.begin(), files.end());

	for(unsigned int i = 0; i< files0.size(); ++i) {
		const std::string &file = files0[i];
		g_debug(" -- init: opening file: " << file << "; " << omp_get_thread_num());
		std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
		liblas::Reader lasReader = rf.CreateWithStream(instr);
		liblas::Header lasHeader = lasReader.GetHeader();
		LASPoint pt;
		while(lasReader.ReadNextPoint()) {
			pt.update(lasReader.GetPoint());
			m_pf->addPoint(pt);
		}
	}

	m_pf->flushAll();
	m_pf->close();
	m_pf->openRead();
}

const Bounds& SortedPointStream::bounds() const {
	return m_bounds;
}

bool SortedPointStream::next(std::list<std::shared_ptr<LASPoint> > &pts) {
	return m_pf->nextRow(pts);
}
