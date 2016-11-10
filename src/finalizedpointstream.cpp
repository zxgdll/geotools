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

	/*
	liblas::ReaderFactory rf;
	Bounds bounds;
	for(const std::string &file : m_files) {
		std::ifstream str(file.c_str(), std::ios::binary);
		liblas::Reader r = rf.CreateWithStream(str);
		liblas::Header h = r.GetHeader();
		bounds.collapse();
		if(!LasUtil::computeLasBounds(h, bounds))
			LasUtil::computeLasBounds(r, bounds);
		m_bounds.extend(bounds);
	}
	m_bounds.snap(m_cellSize);
	g_debug(" -- bounds " << m_bounds.print());

	m_cells.reset(new MemRaster<uint32_t>(cols(), rows(), true));
	m_cells->fill(0);

	for(const std::string &file : m_files) {
		g_debug(" -- reading file " << file);
		std::ifstream str(file.c_str(), std::ios::binary);
		liblas::Reader r(str);
		liblas::Header h = r.GetHeader();
		LASPoint pt;
		while(r.ReadNextPoint()) {
			pt.update(r.GetPoint());
			size_t idx = toIdx(pt);
			m_cells->set(idx, m_cells->get(idx) + 1);
			++m_pointCount;
		}
	}
	*/

	for(const std::string &file : m_files) {
		LASReader r(file);
		m_bounds.extend(r.bounds());
		m_pointCount += r.pointCount();
	}

	m_cells.reset(new MemRaster<uint32_t>(cols(), rows(), true));
	m_cells->fill(0);

	for(const std::string &file : m_files) {
		LASReader r(file);
		LASPoint pt;
		while(r.next(pt)) {
			size_t idx = toIdx(pt);
			m_cells->set(idx, m_cells->get(idx) + 1);
		}
	}
	m_bounds.snap(m_cellSize);
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
	/*
	if(m_fileIdx >= m_files.size())
		return false;
	if(!m_reader.get()) {
		g_debug(" -- reading file " << m_files[m_fileIdx]);
		m_instr.reset(new std::ifstream(m_files[m_fileIdx].c_str(), std::ios::binary));
		m_reader.reset(new liblas::Reader(*m_instr));
	}
	if(!m_reader->ReadNextPoint()) {
		delete m_reader.release();
		delete m_instr.release();
		if(++m_fileIdx >= m_files.size())
			return false;
		g_debug(" -- reading file " << m_files[m_fileIdx]);
		m_instr.reset(new std::ifstream(m_files[m_fileIdx].c_str(), std::ios::binary));
		m_reader.reset(new liblas::Reader(*m_instr));
		if(!m_reader->ReadNextPoint())
			return false;
	}
	pt.update(m_reader->GetPoint());
	size_t idx = toIdx(pt);
	uint32_t count = m_cells->get(idx);
	if(count == 1) {
		m_cells->set(idx, 0);
		*finalIdx = idx;
	} else if(count > 1) {
		m_cells->set(idx, count - 1);
		*finalIdx = 0;
	} else {
		g_runerr("Illegal count.");
	}
	return true;
	*/
	if(m_fileIdx >= m_files.size())
		return false;
	if(!m_reader.get()) {
		g_debug(" -- reading file " << m_files[m_fileIdx]);
		m_reader.reset(new LASReader(m_files[m_fileIdx]));
	}
	if(!m_reader->next(pt)) {
		m_reader.release();
		if(++m_fileIdx >= m_files.size())
			return false;
		m_reader.reset(new LASReader(m_files[m_fileIdx]));
		if(!m_reader->next(pt))
			return false;
	}
	size_t idx = toIdx(pt);
	uint32_t count = m_cells->get(idx);
	if(count == 1) {
		m_cells->set(idx, 0);
		*finalIdx = idx;
	} else if(count > 1) {
		m_cells->set(idx, count - 1);
		*finalIdx = 0;
	} else {
		g_runerr("Illegal count.");
	}
	return true;
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

