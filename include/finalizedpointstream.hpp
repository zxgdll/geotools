#ifndef __FINALIZEDPOINTSTREAM_HPP__
#define __FINALIZEDPOINTSTREAM_HPP__

#include <vector>
#include <unordered_map>
#include <memory>
#include <fstream>
#include <queue>

#include "liblas/liblas.hpp"

#include "laspoint.hpp"
#include "util.hpp"

namespace geotools {
namespace las {

class FinalizedPointStream {
private:
	std::vector<std::string> m_files;
	std::unordered_map<size_t, size_t> m_cells;
	geotools::util::Bounds m_bounds;
	double m_cellSize;
	size_t m_fileIdx;
	std::unique_ptr<std::ifstream> m_instr;
	std::unique_ptr<liblas::Reader> m_reader;
	size_t m_cols;
	size_t m_pointCount;

	void init();

public:
	FinalizedPointStream(const std::vector<std::string> &files, double cellSize);
	bool next(LASPoint &pt, size_t *finalIdx);
	size_t pointCount() const;
	const Bounds& bounds() const;
	size_t toIdx(const LASPoint &pt) const;
	size_t cols() const;
	size_t rows() const;
	size_t toCol(const LASPoint &pt) const;
	size_t toRow(const LASPoint &pt) const;
};

}
}
#endif

