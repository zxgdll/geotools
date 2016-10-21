#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <memory>
#include <cstdio>

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

SortedPointStream::SortedPointStream(const std::list<std::string> &files, double blockSize) {
	init(files, blockSize);
}

SortedPointStream::~SortedPointStream() {
	for(const unsigned int &row : m_rows) {
		try {
			std::remove(_rowFile(row).c_str());
		} catch(...) {
			g_warn("Failed to delete " << _rowFile(row));
		}
	}
}

unsigned int SortedPointStream::pointCount() const {
	return m_pointCount;
}

unsigned int SortedPointStream::rowCount() const {
	return m_rows.size();
}

void SortedPointStream::init(const std::list<std::string> &files, double blockSize) {

	m_bounds.collapse();
	m_pointCount = 0;
	m_row = 0;

	std::map<unsigned int, std::string> rowFileNames;

	for(const std::string &file : files) {
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
	m_bounds.snap(g_abs(blockSize));

	std::map<unsigned int, unsigned int> rowCounts;
	std::map<unsigned int, std::unique_ptr<std::ofstream> > rowFiles;

	for(const std::string &file : files) {
		g_debug(" -- init: opening file: " << file);
		std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
		liblas::ReaderFactory rf;
		liblas::Reader lasReader = rf.CreateWithStream(instr);
		liblas::Header lasHeader = lasReader.GetHeader();

		g_debug(" -- init writing points");
		LASPoint pt;
		while(lasReader.ReadNextPoint()) {
			pt.update(lasReader.GetPoint());
			unsigned int row = _getRow(pt.y, m_bounds, blockSize);
			if(rowFiles.find(row) == rowFiles.end()) {
				std::unique_ptr<std::ofstream> os(new std::ofstream(_rowFile(row), std::ios::binary));
				os->seekp(sizeof(unsigned int), std::ios::beg);
				rowFiles[row] = std::move(os);
				rowCounts[row] = 0;
				m_rows.push_back(row);
			}

			pt.write(*(rowFiles[row].get()));
			rowCounts[row]++;
			++m_pointCount;
		}
	}

	for(const auto &it : rowFiles) {
		unsigned int c = rowCounts[it.first];
		it.second->seekp(0, std::ios::beg);
		it.second->write((char *) &c, sizeof(unsigned int));
	}

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

void SortedPointStream::next(std::list<LASPoint> &pts) {
	#pragma omp critical(__sps_init)
	{
		if(m_row < m_rows.size()) {

			unsigned int row = m_rows[m_row];

			std::ifstream instr(_rowFile(row).c_str(), std::ios::binary);
			unsigned int count;
			instr.read((char *) &count, sizeof(unsigned int));
			if(count) {
				LASPoint pt;
				for(unsigned int i = 0; i < count; ++i) {
					pt.read(instr);
					pts.push_back(pt); // copy
				}
				/*
				double x1 = G_DBL_MAX_POS, y1 = G_DBL_MAX_POS;
				double x2 = G_DBL_MAX_NEG, y2 = G_DBL_MAX_NEG;
				while(count--) {
					pt.read(instr);
					x1 = g_min(pt.x, x1);
					y1 = g_min(pt.y, y1);
					x2 = g_max(pt.x, x2);
					y2 = g_max(pt.y, y2);
					pts.push_back(LASPoint(pt));
				}
				m_boundsTracker.add(x1, y1, x2, y2);
				*/
			}

			++m_row;
		}
	}
}
