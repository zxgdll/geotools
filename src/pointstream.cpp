#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <memory>
#include <cstdio>

#include "liblas/liblas.hpp"

#include "geotools.h"
#include "pointstream.hpp"
#include "lasutil.hpp"

using namespace geotools::las;
using namespace geotools::util;

double LASPoint::scaleX = 0;
double LASPoint::scaleZ = 0;
double LASPoint::scaleY = 0;

PointStream::PointStream(const std::list<std::string> &files, bool deepBounds) {
	init(files, deepBounds);
}

PointStream::~PointStream() {
	delete m_lasReader;
	delete m_instr;
}

unsigned int PointStream::pointCount() {
	return m_pointCount;
}

void PointStream::init(const std::list<std::string> &files, bool deepBounds) {

	m_files.clear();
	m_files.assign(files.begin(), files.end());
	m_fileIdx = 0;

	liblas::ReaderFactory rf;
	m_pointCount = 0;
	m_lasReader = nullptr;
	m_currentFile = "";
	m_bounds.collapse();

	for(const std::string &file : files) {
		g_debug(" -- init: opening file: " << file);
		
		std::ifstream istr(file.c_str(), std::ios::binary);
		liblas::Reader lasReader = rf.CreateWithStream(istr);
		liblas::Header lasHeader = lasReader.GetHeader();

		LASPoint::setScale(lasHeader.GetScaleX(), lasHeader.GetScaleY(), lasHeader.GetScaleZ());
		m_pointCount += lasHeader.GetPointRecordsCount();

		g_debug(" -- init computing bounds");
		Bounds bounds;
		bounds.collapse();
		if(!LasUtil::computeLasBounds(lasHeader, bounds, 2) && deepBounds)
			LasUtil::computeLasBounds(lasReader, bounds, 2); // If the header bounds are bogus.
		m_bounds.extend(bounds);
		m_fileBounds[file] = bounds;
	}

}

Bounds PointStream::bounds() {
	return m_bounds;
}

bool PointStream::contains(double x, double y) {
	for(const auto &it : m_fileBounds) {
		if(it.second.contains(x, y))
			return true;
	}
	return false;
}

bool PointStream::contains(double x1, double y1, double x2, double y2) {
	for(const auto &it : m_fileBounds) {
		if(it.second.contains(x1, y1) && it.second.contains(x2, y2))
			return true;
	}
	return false;
}

bool PointStream::containsCompleted(double x, double y) {
	for(const Bounds &bounds : m_completedBounds) {
		if(bounds.contains(x, y))
			return true;
	}
	return false;
}

bool PointStream::containsCompleted(double x1, double y1, double x2, double y2) {
	for(const Bounds &bounds : m_completedBounds) {
		if(bounds.contains(x1, y1) && bounds.contains(x2, y2))
			return true;
	}
	return false;
}

bool PointStream::loadNextReader() {
	if(m_fileIdx > 0) {
		m_completedBounds.push_back(m_fileBounds[m_currentFile]);
		if(m_fileIdx >= m_files.size())
			return false;
	}
	liblas::ReaderFactory rf;
	m_currentFile = m_files[m_fileIdx];
	m_instr = new std::ifstream(m_currentFile, std::ios::binary);
	m_lasReader = new liblas::Reader(*m_instr);
	m_curPt = 0;
	m_curPtCount = m_lasReader->GetHeader().GetPointRecordsCount();
	++m_fileIdx;
	return true;
}

bool PointStream::next(LASPoint &pt, bool *lastPoint) {
	while(!m_lasReader || !m_lasReader->ReadNextPoint()) {
		if(!loadNextReader()) {
			*lastPoint = true;
			return false;
		}
	}
	pt.update(m_lasReader->GetPoint());
	*lastPoint = ++m_curPt == m_curPtCount;
	return true;
}
