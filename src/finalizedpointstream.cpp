#include <fstream>
#include <thread>
#include <mutex>

#include "liblas/liblas.hpp"

#include "geotools.h"
#include "lasutil.hpp"
#include "finalizedpointstream.hpp"

using namespace geotools::las;

FinalizedPointStream::FinalizedPointStream(const std::vector<std::string> &files, 
		double cellSize, uint32_t threads) :
	m_files(files), 
	m_cellSize(cellSize),
	m_threads(threads) {

	init();
}

void _initializer(FinalizedPointStream *fps) {
	fps->initializer();
}

void FinalizedPointStream::initializer() {
	std::string file;
	while(true) {
		if(!m_fileq.size())
			return;
		m_fqm.lock();
		if(!m_fileq.size()) {
			m_fqm.unlock();
			return;
		}
		file = m_fileq.front();
		m_fileq.pop();
		m_fqm.unlock();

		g_debug(" -- reading file " << file << "; " << m_fileq.size());
		std::ifstream str(file.c_str(), std::ios::binary);
		liblas::Reader r(str);
		liblas::Header h = r.GetHeader();
		LASPoint pt;
		while(r.ReadNextPoint()) {
			pt.update(r.GetPoint());
			size_t idx = toIdx(pt);
			m_fqm.lock();
			if(m_cells.find(idx) == m_cells.end())
				m_cells[idx] = 0;
			++m_cells[idx];
			++m_pointCount;
			m_fqm.unlock();
		}
	}
}

void FinalizedPointStream::init() {

	g_debug(" -- init " << m_cellSize);
	m_bounds.collapse();
	m_pointCount = 0;
	m_fileIdx = 0;

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

	for(const std::string &file : m_files)
		m_fileq.push(file);

	std::list<std::thread> threads;
	g_debug(" -- starting initialization")
	for(uint32_t i = 0; i < m_threads; ++i) {
		std::thread t(_initializer, this);
		threads.push_back(std::move(t));
	}
	for(std::thread &t : threads)
		t.join();
	g_debug(" -- finished initialization")
} 

size_t FinalizedPointStream::pointCount() const {
	return m_pointCount;
}

const Bounds& FinalizedPointStream::bounds() const {
	return m_bounds;
}

bool FinalizedPointStream::next(LASPoint &pt, size_t *finalIdx) {
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
	if(--m_cells[idx] == 0) {
		m_cells.erase(idx);
		*finalIdx = idx;
	} else {
		*finalIdx = 0;
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
	return toRow(pt) * cols() + toCol(pt);
}

