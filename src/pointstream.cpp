#include <fstream>
#include <sstream>
#include <map>
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

PointStream::PointStream(const std::string &file, bool deepBounds) :
	m_file(file) {
	init(deepBounds);
}

PointStream::~PointStream() {
	delete m_lasReader;
	delete m_instr;
}

unsigned int PointStream::pointCount() {
	return m_pointCount;
}

void PointStream::init(bool deepBounds) {
	g_debug(" -- init: opening file: " << m_file);
	
	m_instr = new std::ifstream (m_file.c_str(), std::ios::binary);

	liblas::ReaderFactory rf;
	liblas::Reader lasReader = rf.CreateWithStream(*m_instr);
	liblas::Header lasHeader = lasReader.GetHeader();

	LASPoint::setScale(lasHeader.GetScaleX(), lasHeader.GetScaleY(), lasHeader.GetScaleZ());
	m_pointCount = lasHeader.GetPointRecordsCount();

	g_debug(" -- init computing bounds");
	if(!LasUtil::computeLasBounds(lasHeader, m_fileBounds, 2) && deepBounds)
		LasUtil::computeLasBounds(lasReader, m_fileBounds, 2); // If the header bounds are bogus.

	lasReader.Reset();

	m_lasReader = new liblas::Reader(*m_instr);
}

Bounds PointStream::fileBounds() {
	return m_fileBounds;
}

bool PointStream::next(LASPoint &pt) {
	if(m_lasReader->ReadNextPoint()) {
		pt.update(m_lasReader->GetPoint());
		return true;
	}
	return false;
}
