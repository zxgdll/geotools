#include <fstream>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "liblas/liblas.hpp"

#include "geotools.h"
#include "lasutil.hpp"
#include "finalizedpointstream.hpp"
#include "raster.hpp"

using namespace geotools::las;
using namespace geotools::raster;

FinalizedPointStream::FinalizedPointStream(const std::vector<std::string> &files, 
		double cellSize) :
	m_files(files), 
	m_cellSize(cellSize) {

	init();
}

void FinalizedPointStream::init() {

	g_debug(" -- init " << m_cellSize);
	m_bounds.collapse();
	m_pointCount = 0;
	m_fileIdx = 0;

	for(const std::string &file : m_files) {
		LASReader r(file);
		m_bounds.extend(r.bounds());
		m_pointCount += r.pointCount();
	}
	m_bounds.snap(m_cellSize);

	m_cells.reset(new MemRaster<uint32_t>(cols(), rows(), true));
	m_cells->fill(0);

	for(const std::string &file : m_files) {
		g_debug(" -- counting " << file);
		LASReader r(file);
		LASPoint pt;
		while(r.next(pt)) {
			size_t idx = toIdx(pt);
			m_cells->set(idx, m_cells->get(idx) + 1);
		}
	}
	g_debug(" -- bounds " << m_bounds.print());
	g_debug(" -- finished initialization")
} 

size_t FinalizedPointStream::pointCount() const {
	return m_pointCount;
}

const Bounds& FinalizedPointStream::bounds() const {
	return m_bounds;
}

bool FinalizedPointStream::next(LASPoint &pt, size_t *finalIdx) {
	if(!m_reader.get()) {
		g_debug(" -- reading file 1 " << m_files[m_fileIdx]);
		m_reader.reset(new LASReader(m_files[m_fileIdx]));
		++m_fileIdx;
	}
	if(!m_reader->next(pt)) {
		delete m_reader.release();
		if(m_fileIdx >= m_files.size())
			return false;
		g_debug(" -- reading file 2 " << m_files[m_fileIdx]);
		m_reader.reset(new LASReader(m_files[m_fileIdx]));
		++m_fileIdx;
		if(!m_reader->next(pt))
			return false;
	}
	size_t idx = toIdx(pt);
	uint32_t count = m_cells->get(idx);
	if(count > 0) {
		m_cells->set(idx, count - 1);
		*finalIdx = count == 1 ? idx : 0;
		return true;
	}
	return false;
}

size_t FinalizedPointStream::cols() const {
	return (size_t) (m_bounds.width() / m_cellSize) + 1;
}

size_t FinalizedPointStream::rows() const {
	return (size_t) (m_bounds.height() / m_cellSize) + 1;
}

size_t FinalizedPointStream::toCol(const LASPoint &pt) const {
	return (size_t) ((pt.x - m_bounds.minx()) / m_cellSize);
}

size_t FinalizedPointStream::toRow(const LASPoint &pt) const {
	return rows() - (size_t) ((pt.y - m_bounds.miny()) / m_cellSize) - 1;
}

size_t FinalizedPointStream::toIdx(const LASPoint &pt) const {
	//_debug(" -- toidx " << pt.x << ", " << pt.y << "; " << toRow(pt) << "; " << cols() << "; " << toCol(pt));
	return toRow(pt) * cols() + toCol(pt);
}

