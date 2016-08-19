/*
 * Raster.hpp
 *
 *  Created on: Jan 21, 2016
 *      Author: rob
 */

#ifndef INCLUDE_RASTER_HPP_
#define INCLUDE_RASTER_HPP_

#include <queue>
#include <stdexcept>
#include <map>
#include <unordered_map>
#include <vector>
#include <cstring>

#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <eigen3/Eigen/Core>

#include "geotools.h"
#include "util.hpp"

using namespace geotools::util;

namespace geotools {

namespace raster {

/**
 * Simple class to represent a single grid cell.
*/
class Cell {
public:
	int col;
	int row;
	Cell(int col, int row) :
		col(col), row(row) {
	}
};

/**
 * Used by Grid::floodFill to determine whether
 * a pixel should be filled.
 */
template <class T>
class FillOperator {
public:
	virtual bool fill(T value) const =0;
};

template <class T>
class TargetOperator : public FillOperator<T> {
private:
	T m_match;
public:
	TargetOperator(T match) :
		m_match(match) {
	}

	bool fill(T value) const {
		return value == m_match;
	}
};

/**
 * Abstract class for grids (rasters).
 */
template <class T>
class Grid {
protected:
	T m_min = 0;
	T m_max = 0;
	T m_mean = 0;
	T m_stddev = 0;
	T m_variance = 0;
	T m_sum = 0;
	int m_count = 0;
	bool m_stats = false;

	/**
	 * Compute the table of Gaussian weights given the size of the table
	 * and the std. deviation.
	 */
	void gaussianWeights(double *weights, int size, double sigma) {
		if(size % 2 == 0) ++size;
		for(int r = 0; r < size; ++r) {
			for(int c = 0; c < size; ++c) {
				int x = size / 2 - c;
				int y = size / 2 - r;
				weights[r * size + c] = (1 / (2 * G_PI * sigma * sigma)) * pow(G_E, -((x * x + y * y) / (2.0 * sigma * sigma)));
			}
		}
	}

public:
	
	virtual int rows() const =0;

	virtual int cols() const =0;
	
	virtual size_t size() const =0;

	virtual void fill(const T value) =0;

	/**
	 * Return a pointer to an in-memory grid of the raster data.
	 * Throw an exception if this is not possible.
	 */
	virtual T *grid() =0;
	
	/**
	 * Returns true if this class has a complete, in-memory
	 * grid that can be manipulated.
	 */
	virtual bool hasGrid() const =0;

	/**
	 * Return a reference to the value held at
	 * the given index in the grid.
	 * Not const because the get operation might imply (e.g.)
	 * a buffering operation in the subclass.
	 */
	virtual T &get(size_t idx) =0;

	virtual T &get(int col, int row) =0;

	virtual void set(int col, int row, const T value) =0;

	virtual void set(size_t idx, const T value) =0;

	virtual bool has(int col, int row) const =0;
	virtual bool has(size_t idx) const =0;

	virtual bool isNoData(int col, int row) =0;
	virtual bool isNoData(size_t idx) =0;

	/**
	 * Return the element at the given index.
	 */
	virtual T &operator[](size_t idx) =0;

	virtual bool isSquare() const =0;

	virtual T nodata() const =0;

	virtual void nodata(T nodata) =0;

	/**
	 * Load the contents of a single block into the given Grid instance.
	 */
	virtual void readBlock(int col, int row, Grid<T> &block, int srcCol = 0, int srcRow = 0) =0;

	/**
	 * Write the given block to the raster.
	 */
	virtual void writeBlock(int col, int row, Grid<T> &block, int dstCol = 0, int dstRow = 0) =0;

	/**
	 * Write the full block to the raster.
	 */
	virtual void writeBlock(Grid<T> &block) =0;

	/**
	 * Read the full block from the raster.
	 */
	virtual void readBlock(Grid<T> &block) =0;

	void computeStats() {
		size_t i;
		for(i = 0; i < size(); ++i) {
			if(!isNoData(i)) {
				m_min = m_max = get(i);
				break;
			}
		}

		m_sum = 0;
		m_count = 0;
		T m = 0;
		T s = 0;
		int k = 1;
		// Welford's method for variance.
		// i has the index of the first non-nodata element.
		for(; i < size(); ++i) {
			T v = get(i);
			if(v != nodata()) {
				T oldm = m;
				m = m + (v - m) / k;
				s = s + (v - m) * (v - oldm);
				m_sum += v;
				m_min = g_min(m_min, v);
				m_max = g_max(m_max, v);
				++m_count;
				++k;
			}
		}
		m_mean = m_sum / m_count;
		m_variance = s / m_count;
		m_stddev = std::sqrt(m_variance);
		m_stats = true;
		g_trace("Count: " << m_count << "; Sum: " << m_sum << "; Min: " << m_min 
			<< "; Max: " << m_max << "; Mean: " << m_mean 
			<< "; Variance: " << m_variance << "; Std Dev: " << m_stddev);
	}

	T max() {
		if(!m_stats)
			computeStats();
		return m_max;
	}

	T min() {
		if(!m_stats)
			computeStats();
		return m_min;
	}

	T mean() {
		if(!m_stats)
			computeStats();
		return m_mean;
	}
	
	T stddev() {
		if(!m_stats)
			computeStats();
		return m_stddev;
	}

	T variance() {
		if(!m_stats)
			computeStats();
		return m_variance;
	}

	/**
	 * Fill the grid, beginning with the target cell, where any contiguous cell
	 * satisfies the given FillOperator. The other grid is actually filled,
	 * and the present grid is unchanged *unless* the present grid is passed
	 * as other.
	 */
	template <class U>
	std::vector<int> floodFill(int col, int row, FillOperator<T> &op, Grid<U> &other, U fill) {

		int minc = cols() + 1;
		int minr = rows() + 1;
		int maxc = -1;
		int maxr = -1;
		int area = 0;

		std::queue<Cell*> q;
		q.push(new Cell(col, row));

		std::vector<bool> visited(size(), false); // Tracks visited pixels.

		while(q.size()) {
			
			Cell *cel = q.front();
			row = cel->row;
			col = cel->col;
			q.pop();
			delete cel;

			size_t idx = (size_t) row * cols() + col;
			
			if(!visited[idx] && op.fill(get(col, row))) {

				minc = g_min(col, minc);
				maxc = g_max(col, maxc);
				minr = g_min(row, minr);
				maxr = g_max(row, maxr);
				++area;
				other.set(col, row, fill);
				visited[idx] = true;

				if(row > 0)
					q.push(new Cell(col, row - 1));
				if(row < rows() - 1)
					q.push(new Cell(col, row + 1));
				
				for(int c = col - 1; c >= 0; --c) {
					idx = (size_t) row * cols() + c;
					if(!visited[idx] && op.fill(get(c, row))) {
						minc = g_min(c, minc);
						++area;
						other.set(c, row, fill);
						visited[idx] = true;
						if(row > 0)
							q.push(new Cell(c, row - 1));
						if(row < rows() - 1)
							q.push(new Cell(c, row + 1));
					} else {
						break;
					}
				}	
				for(int c = col + 1; c < cols(); ++c) {
					idx = (size_t) row * cols() + c;
					if(!visited[idx] && op.fill(get(c, row))) {
						maxc = g_max(c, maxc);
						++area;
						other.set(c, row, fill);
						visited[idx] = true;
						if(row > 0)
							q.push(new Cell(c, row - 1));
						if(row < rows() - 1) 
							q.push(new Cell(c, row + 1));
					} else {
						break;
					}
				}	
			}
		}

		return {minc, minr, maxc, maxr, area};
	}

	std::vector<int> floodFill(int col, int row, T target, T fill) {
		TargetOperator<T> op(target);
		return floodFill(col, row, op, *this, fill);
	}

	std::vector<int> floodFill(int col, int row, FillOperator<T> &op, T fill) {
		return floodFill(col, row, op, *this, fill);
	}


	/**
	 * Smooth the raster and write the smoothed version to the output raster.
	 */
	// TODO: No accounting for nodata.
	void smooth(Grid<float> &smoothed, double sigma, int size) {
		if(sigma <= 0)
			g_argerr("Sigma must be > 0.");
		if(size < 3)
			g_argerr("Kernel size must be 3 or larger.");
		if(size % 2 == 0) {
			g_warn("Kernel size must be odd. Rounding up.");
			size++;
		}
		double weights[size * size];
		gaussianWeights(weights, size, sigma);
		for(int r = 0; r < rows() - size; ++r) {
			for(int c = 0; c < cols() - size; ++c) {
				double t = 0.0;
				double v;
				bool foundNodata = false;
				for(int gr = 0; gr < size; ++gr) {
					for(int gc = 0; gc < size; ++gc) {
						v = get(c + gc, r + gr);
						if(v == nodata()) {
							foundNodata = true;
							break;
						}
						t += weights[gr * size + gc] * v;
					}
					if(foundNodata)
						break;
				}
				if(!foundNodata)
					smoothed.set(c + size / 2, r + size / 2, t);
			}
		}
	}
	
	/**
	 * The radius is given with cells as the unit, but
	 * can be rational. When determining which cells to
	 * include in the calculation, any cell which partially
	 * falls in the radius will be included.
	 */
	void voidFillIDW(double radius, int count = 4, double exp = 2.0);

};

/**
 * A convenience class for managing a grid of values.
 * Handles allocation and deallocation of memory.
 */
template <class T>
class MemRaster : public Grid<T> {
private:
	int m_cols;
	int m_rows;
	T *m_grid;
	T m_nodata;

	void (*m_item_dealloc)(T);

	/**
	 * Checks if the grid has been initialized. Throws exception otherwise.
	 */
	void checkInit() const {
		if(m_grid == nullptr)
			g_runerr("This instance has not been initialized.");
	}

public:
	MemRaster() :
		m_grid(nullptr),
		m_cols(-1), m_rows(-1),
		m_item_dealloc(nullptr),
		m_nodata(0) {
	}

	MemRaster(int cols, int rows) : MemRaster() {
		init(cols, rows);
	}

	template <class U>
	MemRaster(Grid<U> &tpl) : MemRaster() {
		init(tpl.cols(), tpl.rows());
	}

	~MemRaster() {
		if(m_item_dealloc != nullptr) {
			for(size_t i = 0; i < (size_t) m_cols * m_rows; ++i)
				m_item_dealloc(m_grid[i]);
		}
		if(m_grid != nullptr)
			free(m_grid);
	}

	/**
	 * A pointer to a function that can deallocate the grid
	 * items. Not necessary if a primitive type is used (the usual case.)
	 */
	void setDeallocator(void (*item_dealloc)(T)) {
		m_item_dealloc = item_dealloc;
	}

	/**
	 * Return a pointer to the allocated memory.
	 */
	T *grid() {
		return m_grid;
	}

	bool hasGrid() const {
		return true;
	}

	/**
	 * Cas a MemRaster to some other type.
	 */
	template <class U>
	void convert(MemRaster<U> &g) {
		g.init(cols(), rows());
		for(size_t i = 0; i < size(); ++i)
			g.set(i, (U) get(i));
	}

	int rows() const {
		return m_rows;
	}

	int cols() const {
		return m_cols;
	}

	size_t size() const {
		return (size_t) m_rows * m_cols;
	}

	template <class U>
	void init(Grid<U> &tpl) {
		init(tpl.cols(), tpl.rows());
	}

	/**
	 * Initialize with the given number of cols and rows.
	 * (Re)allocates memory for the internal grid.
	 */
	void init(int cols, int rows) {
		if(cols <= 0 || rows <= 0)
			g_argerr("Invalid row or column count.");
		if(cols != m_cols || rows != m_rows) {
			m_cols = cols;
			m_rows = rows;
			if(m_grid != nullptr)
				free(m_grid);
			m_grid = (T *) malloc(sizeof(T) * cols * rows);
		}
	}

	/**
	 * Fill the grid with the given value.
	 */
	void fill(const T value) {
		checkInit();
		for(size_t i = 0; i < size(); ++i)
			m_grid[i] = value;	
	}

	/**
	 * Return a reference to the value held at
	 * the given index in the grid.
	 */
	T &get(size_t idx) {
		checkInit();
		if(idx >= size())
			g_argerr("Index out of bounds.");
		return m_grid[idx];
	}

	T &get(int col, int row) {
		size_t idx = (size_t) row * m_cols + col;
		return get(idx);
	}

	bool isNoData(int col, int row) {
		return get(col, row) == m_nodata;
	}

	bool isNoData(size_t idx) {
		return get(idx) == m_nodata;
	}

	void set(int col, int row, const T value) {
		size_t idx = (size_t) row * m_cols + col;
		set(idx, value);
	}

	void set(size_t idx, const T value) {
		checkInit();
		if(idx >= size())
			g_argerr("Index out of bounds.");
		m_grid[idx] = value;
	}

	bool has(int col, int row) const {
		return col >= 0 && col < m_cols && row >= 0 && row < m_rows;
	}

	bool has(size_t idx) const {
		return idx < (size_t) m_cols * m_rows;
	}

	/**
	 * Return the element at the given index.
	 */
	T &operator[](size_t idx) {
		checkInit();
		if(idx >= size())
			g_argerr("Index out of bounds.");
		return m_grid[idx];
	}

	bool isSquare() const {
		return cols() == rows();
	}

	void toMatrix(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx) {
		for(int r = 1; r < rows(); ++r) {
			for(int c = 0; c < cols(); ++c)
				mtx(r, c) = get(c, r);
		}
	}

	void fromMatrix(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx) {
		for(int r = 1; r < rows(); ++r) {
			for(int c = 0; c < cols(); ++c)
				set(c, r, mtx(r, c));
		}
	}

	T nodata() const {
		return m_nodata;
	}

	void nodata(T nodata) {
		m_nodata = nodata;
	}

	/**
	 * Load the contents of a single block into the given Grid instance.
	 */
	void readBlock(int col, int row, Grid<T> &block, int dstCol = 0, int dstRow = 0) {
		if(&block == this)
			g_argerr("Recursive call to readBlock.");
		if(col + block.cols() - dstRow > m_cols)
			g_argerr("Block is wider than the available space.");
		if(row + block.rows() - dstCol > m_rows)
			g_argerr("Block is taller than the available space.");
		if(block.hasGrid()) {
			for(int r = 0; r < block.rows() - dstRow; ++r) {
				std::memcpy(
					block.grid() + (r + dstRow) * block.cols() + dstCol, 
					m_grid + (row + r) * m_cols + col, 
					(block.cols() - dstCol) * sizeof(T)
				);
			}
		} else {
			for(int r = 0; r < block.rows() - dstRow; ++r) {
				for(int c = 0; c < block.cols() - dstCol; ++c)
					block.set(c + dstCol, r + dstRow, get(c + col, r + row));
			}
		}
	}

	/**
	 * Write a part of the given block to the raster.
	 */
	void writeBlock(int col, int row, Grid<T> &block, int srcCol = 0, int srcRow = 0) {
		if(&block == this)
			g_argerr("Recursive call to writeBlock.");
		if(col + block.cols() - srcCol > m_cols)
			g_argerr("Block is wider than the available space.");
		if(row + block.rows() - srcRow > m_rows)
			g_argerr("Block is taller than the available space.");
		if(block.hasGrid()) {
			for(int r = 0; r < block.rows(); ++r) {
				std::memcpy(
					m_grid + (row + r) * m_cols + col, 
					block.grid() + (r + srcRow) * block.cols() + srcCol, 
					(block.cols() - srcRow) * sizeof(T)	
				);
			}
		} else {
			for(int r = 0; r < block.rows() - srcRow; ++r) {
				for(int c = 0; c < block.cols() - srcCol; ++c)
					set(c + col, r + row, block.get(c + srcCol, r + srcRow));
			}
		}
	}

	/**
	 * Write the full block to the raster.
	 */
	void writeBlock(Grid<T> &block) {
		writeBlock(0, 0, block);
	}

	/**
	 * Read the full block from the raster.
	 */
	void readBlock(Grid<T> &block) {
		readBlock(0, 0, block);
	}


};

template <class T>
class BlockCache {
private: 
	size_t m_size;
	int m_bw;
	int m_bh;
	GDALRasterBand *m_band;
	std::unordered_map<size_t, T*> m_blocks;       // idx, block
	std::unordered_map<size_t, size_t> m_idx_time; // idx, time
	std::unordered_map<size_t, bool> m_dirty;      // idx, bool
	std::map<size_t, size_t> m_time_idx;           // time, idx
	size_t m_time;

	void flushBlock(size_t idx) {
		if(hasBlock(idx) && m_band->GetDataset()->GetAccess() == GA_Update) {
			T *blk = m_blocks[idx];
			if(m_band->WriteBlock(toCol(idx) / m_bw, toRow(idx) / m_bh, blk) != CE_None)
				g_runerr("Failed to flush block.");
			m_dirty[idx] = false;
		}
	}

	size_t toIdx(int col, int row) {
		//g_trace("to idx: " << col << "," << row << " --> " << (((size_t) (col / m_bw) << 32) | (row / m_bh)));
		return ((size_t) (col / m_bw) << 32) | (row / m_bh);
	}

	int toCol(size_t idx) {
		//g_trace("to col: " << idx << " --> " << ((idx >> 32) & 0xffffffff) * m_bw);
		return ((idx >> 32) & 0xffffffff) * m_bw;
	}

	int toRow(size_t idx) {
		//g_trace("to row: " << idx << " --> " << (idx & 0xffffffff) * m_bh);
		return (idx & 0xffffffff) * m_bh;
	}

 	T* freeOldest() {
		auto it = m_time_idx.rbegin();
		size_t time = it->first;
		size_t idx = it->second;
		if(m_dirty[idx])
			flushBlock(idx);
		T *blk = m_blocks[idx];
		m_blocks.erase(idx);
		m_idx_time.erase(idx);
		m_time_idx.erase(time);
		m_dirty.erase(idx);
		return blk;
	}

	T* freeOne() {
		T *blk = nullptr;
		while(m_blocks.size() >= m_size) {
			if(blk)
				free(blk);
			blk = freeOldest();
		}
		return blk;
	}

public:
	BlockCache() :
		m_band(nullptr), 
		m_size(0),
		m_time(0),
		m_bw(0), m_bh(0) {
	}

	int blockWidth() {
		return m_bw;
	}

	int blockHeight() {
		return m_bh;
	}

	size_t toBlockIdx(int col, int row) {
		return (row % m_bh) * m_bw + (col % m_bw);
	}

	void setRasterBand(GDALRasterBand *band) {
		m_band = band;
		band->GetBlockSize(&m_bw, &m_bh);
	}

	bool hasBlock(int col, int row) {
		return hasBlock(toIdx(col, row));
	}

	bool hasBlock(size_t idx) {
		return m_blocks.find(idx) != m_blocks.end();
	}

	void setSize(size_t size) {
		while(m_blocks.size() > size)
			freeOne();
		m_size = size;
	}

	size_t getSize() {
		return m_size;
	}

	T* getBlock(int col, int row, bool forWrite) {
		size_t idx = toIdx(col, row);
		if(!hasBlock(idx)) {
			T *blk = freeOne();
			if(!blk)
				blk = (T *) malloc(sizeof(T) * m_bw * m_bh);
			if(m_band->ReadBlock(col / m_bw, row / m_bh, blk) != CE_None)
				g_runerr("Failed to read block.");
			m_blocks[idx] = blk;
		}
		++m_time; // TODO: No provision for rollover
		m_time_idx.erase(m_idx_time[idx]);
		m_time_idx[m_time] = idx;
		m_idx_time[idx] = m_time;
		if(forWrite)
			m_dirty[idx] = true;
		return m_blocks[idx];
	}

	void flush() {
		for(auto it = m_blocks.begin(); it != m_blocks.end(); ++it)
			flushBlock(it->first);
	}

	void close() {
		flush();
		for(auto it = m_blocks.begin(); it != m_blocks.end(); ++it) {
			if(it->second)			
				free(it->second);
		}
	}

};

template <class T>
class Raster : public Grid<T> {
private:
	int m_cols, m_rows;		// Raster cols/rows
	int m_bandn;			// The band number
	bool m_writable;		// True if the raster is writable
	T m_nodata;			// Nodata value.
	T *m_block;			// Block storage
	GDALDataset *m_ds;		// GDAL dataset
	GDALRasterBand *m_band;		// GDAL band
	double m_trans[6];		// Raster transform
	bool m_inited = false;		// True if the instance is initialized.
	GDALDataType m_type;		// GDALDataType -- limits the possible template types.
	std::string m_filename;		// Raster filename
	BlockCache<T> m_cache;		// Block cache.

	/**
	 * Loads the block that contains the given row and column.
	 */
	void loadBlock(int col, int row, bool forWrite) {
		if(!m_inited)
			g_runerr("Not inited before attempted read.");
		if(!has(col, row))
			g_argerr("Row or column out of bounds:" << col << ", " << row);
		m_block = m_cache.getBlock(col, row, forWrite);
		if(!m_block)
			g_runerr("Failed to load block from cache.");
	}	

	GDALDataType getType(double v) {
		(void) v;
		return GDT_Float64;
	}

	GDALDataType getType(float v) {
		(void) v;
		return GDT_Float32;
	}

	GDALDataType getType(unsigned int v) {
		(void) v;
		return GDT_UInt32;
	}

	GDALDataType getType(int v) {
		(void) v;
		return GDT_Int32;
	}

	GDALDataType getType(unsigned short v) {
		(void) v;
		return GDT_UInt16;
	}

	GDALDataType getType(short v) {
		(void) v;
		return GDT_Int16;
	}

	GDALDataType getType(char v) {
		(void) v;
		return GDT_Byte;
	}

	GDALDataType getType() {
		return getType((T) 0);
	}

	T getDefaultNodata() {
		switch(getType()) {
		case GDT_Float32:
		case GDT_Float64:
			return -9999.0;
		default:
			return 0;
		}
	}

public:

	/**
	 * Basic constructor.
	 */
	Raster() :
		m_cols(-1), m_rows(-1),
		m_bandn(1),
		m_writable(false), 
		m_ds(nullptr), m_band(nullptr), m_block(nullptr),
		m_type(getType()) {
	}

	/**
	 * Create a new raster for writing with a template.
	 */
	template <class D>
	Raster(const std::string &filename, int band, const Raster<D> &tpl) : Raster() {
		std::string proj;
		tpl.projection(proj);
		init(filename, band, tpl.minx(), tpl.miny(), tpl.maxx(), tpl.maxy(), tpl.resolutionX(), 
			tpl.resolutionY(), (T) tpl.nodata(), proj);
	}

	Raster(const std::string &filename, int band, const Raster<T> &tpl) : Raster() {
		init(filename, band, tpl);
	}

	/**
	 * Build a new raster with the given filename, bounds, resolution, nodata and projection.
	 */
	Raster(const std::string &filename, int band, double minx, double miny, double maxx, double maxy,
			double resolutionX, double resolutionY, double nodata, const std::string &proj) : Raster() {
		init(filename,band,  minx, miny, maxx, maxy, resolutionX, resolutionY, nodata, proj);
	}

	Raster(const std::string &filename, int band, Bounds &bounds, double resolutionX, double resolutionY, double nodata, int crs) :
		Raster(filename, band, bounds.minx(), bounds.miny(), bounds.maxx(), bounds.maxy(), 
			resolutionX, resolutionY, nodata, crs) {
	}
	
	/**
	 * Build a new raster with the given filename, bounds, resolution, nodata and SRID.
	 */
	Raster(const std::string &filename, int band, double minx, double miny, double maxx, double maxy,
			double resolutionX, double resolutionY, double nodata, int crs) : Raster() {
		std::string proj = epsg2ProjText(crs);
		init(filename, band, minx, miny, maxx, maxy, resolutionX, resolutionY, nodata, proj);
	}

	/**
	 * Open the given raster and load the given band. Set the writable argument to true
	 * to enable writing.
	 */
	Raster(const std::string &filename, int band = 1, bool writable = false) : Raster() {
		init(filename, band, writable);
	}

	/**
	 * Initialize the raster using the given file (which may not exist) using
	 * another raster as a template. The raster pixel types need not be the same.
	 */
	template <class D>
	void init(const std::string &filename, int band, const Raster<D> &tpl) {
		g_trace("Raster init: " << filename << "; [tpl]");
		std::string proj;
		tpl.projection(proj);
		init(filename, band, tpl.minx(), tpl.miny(), tpl.maxx(), tpl.maxy(), tpl.resolutionX(), 
			tpl.resolutionY(), tpl.nodata(), proj);
	}

	void init(const std::string &filename, int band, const Bounds &bounds, double resolutionX, double resolutionY,
		double nodata, const std::string &proj) {
		init(filename, band, bounds.minx(), bounds.miny(), bounds.maxx(), bounds.maxy(), 
			resolutionX, resolutionY, nodata, proj);
	}

	/**
	 * Initializes the raster with the given filename, bounds, resolution, nodata and projection.
	 */
	void init(const std::string &filename, int band, double minx, double miny, double maxx, double maxy,
			double resolutionX, double resolutionY, double nodata, const std::string &proj) {
		
		g_trace("Raster init: " << filename << ", " << minx << ", " << miny << ", " << maxx << ", " << maxy << ", " << resolutionX << ", " << resolutionY << ", " << nodata << ", " << proj);

		if(resolutionX == 0 || resolutionY == 0)
			g_argerr("Resolution must be larger or smaller than 0.");
		if(maxx < minx)
			g_argerr("Minimum x must be smaller than or equal to maximum x.");
		if(maxy < miny)
			g_argerr("Minimum y must be smaller than or equal to maximum y.");
		if(filename.empty())
			g_argerr("Filename must be given.");

		m_filename.assign(filename);

		// Compute columns/rows
		int width = (int) ((maxx - minx) / resolutionX) * (resolutionX < 0 ? -1 : 1);
		int height = (int) ((maxy - miny) / resolutionY) * (resolutionY < 0 ? -1 : 1);

		// Create GDAL dataset.
		GDALAllRegister();
		m_ds = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(filename.c_str(),
				width, height, 1, m_type, NULL);
		if(m_ds == nullptr)
			g_runerr("Failed to create file.");

		// Initialize geotransform.
		auto lst = std::initializer_list<double>({ 
			resolutionX < 0 ? maxx : minx, resolutionX, 0.0, 
			resolutionY < 0 ? maxy : miny, 0.0, resolutionY 
		});
		std::copy(lst.begin(), lst.end(), m_trans);
		m_ds->SetGeoTransform(m_trans);

		// Set projection.
		if(!proj.empty())
			m_ds->SetProjection(proj.c_str());

		// Save some dataset properties.
		m_rows = m_ds->GetRasterYSize();
		m_cols = m_ds->GetRasterXSize();
		m_band = m_ds->GetRasterBand(band);
		if(m_band == NULL)
			g_runerr("Failed to get band.");
		m_bandn = band;
		m_band->SetNoDataValue(nodata);
		m_nodata = m_band->GetNoDataValue();
		m_cache.setSize(100);
		m_cache.setRasterBand(m_band);
		m_writable = true;
		m_inited = true;
	}

	/**
	 * Initializes a Raster from the existing file.
	 */
	void init(const std::string &filename, int band = 1, bool writable = false) {

		if(filename.empty())
			g_argerr("Filename must be given.");
		
		m_filename.assign(filename);

		// Attempt to open the dataset.
		GDALAllRegister();
		m_ds = (GDALDataset *) GDALOpen(filename.c_str(), writable ? GA_Update : GA_ReadOnly);
		if(m_ds == NULL)
			g_runerr("Failed to open raster.");

		// Save some raster
		m_bandn = band;
		m_ds->GetGeoTransform(m_trans);
		m_band = m_ds->GetRasterBand(band);
		if(m_band == nullptr)
			g_runerr("Failed to get band.");
		m_rows = m_ds->GetRasterYSize();
		m_cols = m_ds->GetRasterXSize();
		m_nodata = m_band->GetNoDataValue();
		m_cache.setSize(100);
		m_cache.setRasterBand(m_band);
		m_writable = writable;
		m_inited = true;
	}

	void setCacheSize(size_t size) {
		m_cache.setSize(size);
	}

	std::string filename() const {
		return m_filename;
	}

	/**
	 * Converts a numerical (EPSG) crs code to a projection string.
	 */
	std::string epsg2ProjText(int crs) {
		OGRSpatialReference ref;
		char *wkt;
		ref.importFromEPSG(crs);
		ref.exportToWkt(&wkt);
		return std::string(wkt);
	}

	/**
	 * Returns true if the raster is initialized.
	 */
	bool inited() {
		return m_inited;
	}

	void fill(T value) {
		MemRaster<T> grd(m_cache.blockWidth(), m_cache.blockHeight());
		grd.fill(value);
		for(int r = 0; r < rows() / grd.rows(); ++r) {
			for(int c = 0; c < cols() / grd.cols(); ++c) {
				if(m_band->WriteBlock(c, r, grd.grid()) != CE_None)
					g_runerr("Fill error.");
			}
		}
	}

	/**
	 * Load the contents of a single block into the given Grid instance.
	 */
	void readBlock(int col, int row, Grid<T> &grd, int dstCol = 0, int dstRow = 0) {
		if(&grd == this)
			g_runerr("Recursive call to readBlock.");
		int cols = g_min(m_cols - col, grd.cols() - dstCol);
		int rows = g_min(m_rows - row, grd.rows() - dstRow);
		if(cols < 1 || rows < 1)
			g_argerr("Zero read size.");
		MemRaster<T> mr(cols, rows);
		if(m_band->RasterIO(GF_Read, col, row, cols, rows, mr.grid(), cols, rows, getType(), 0, 0) != CE_None)
			g_runerr("Failed to read from band.");
		grd.writeBlock(dstCol, dstRow, mr);
	}

	/**
	 * Write a part of the given block to the raster.
	 */
	void writeBlock(int col, int row, Grid<T> &grd, int srcCol = 0, int srcRow = 0) {
		if(&grd == this)
			g_runerr("Recursive call to writeBlock.");
		int cols = g_min(m_cols - col, grd.cols() - srcCol);
		int rows = g_min(m_rows - row, grd.rows() - srcRow);
		if(cols < 1 || rows < 1)
			g_argerr("Zero write size.");
		MemRaster<T> mr(cols, rows);
		grd.readBlock(srcCol, srcRow, mr);
		if(m_band->RasterIO(GF_Write, col, row, cols, rows, mr.grid(), cols, rows, getType(), 0, 0) != CE_None)
			g_runerr("Failed to write to band.");
	}

	/**
	 * Write the full block to the raster.
	 */
	void writeBlock(Grid<T> &block) {
		writeBlock(0, 0, block);
	}
	
	/**
	 * Read the full block from the raster.
	 */
	void readBlock(Grid<T> &block) {
		readBlock(0, 0, block);
	}
	
	/**
	 * Get the x resolution.
	 */
	double resolutionX() const {
		return m_trans[1];
	}

	/**
	 * Get the y resolution.
	 */
	double resolutionY() const {
		return m_trans[5];
	}

	/**
	 * Write the projection data to the given string object.
	 */
	void projection(std::string &proj) const {
		proj.assign(m_ds->GetProjectionRef());
	}

	/**
	 * Return the GDAL datatype of the raster.
	 */
	GDALDataType type() const {
		return m_band->GetRasterDataType();
	}

	/**
	 * The minimum bounding x.
	 */
	double minx() const {
		return toX(0);
	}

	/**
	 * The maximum bounding x.
	 */
	double maxx() const {
		return toX(cols());
	}

	/**
	 * The minimum bounding y.
	 */
	double miny() const {
		return toY(rows());
	}

	/**
	 * The maximum bounding y.
	 */
	double maxy() const {
		return toY(0);
	}

	/**
	 * The width of the raster in map units.
	 */
	double width() const {
		return maxx() - minx();
	}

	/**
	 * The height of the raster in map units.
	 */
	double height() const {
		return maxy() - miny();
	}

	/**
	 * The nodata value.
	 */
	T nodata() const {
		return m_nodata;
	}

	void nodata(T nodata) {
		m_band->SetNoDataValue(nodata);
		m_nodata = nodata;
	}

	/**
	 * Returns the total number of columns.
	 */
	int cols() const {
		return m_cols;
	}

	/**
	 * Returns the total number of rows.
	 */
	int rows() const {
		return m_rows;
	}

	/**
	 * Returns the row for a given y-coordinate.
	 */
	int toRow(double y) const {
		return (int) ((y - m_trans[3]) / m_trans[5]);
	}

	/**
	 * Returns the column for a given x-coordinate.
	 */
	int toCol(double x) const {
		return (int) ((x - m_trans[0]) / m_trans[1]);
	}

	/**
	 * Returns the x-coordinate for a given column.
	 */
	double toX(int col) const {
		return (col * m_trans[1]) + m_trans[0];
	}

	/**
	 * Returns the y-coordinate for a given row.
	 */
	double toY(int row) const {
		return (row * m_trans[5]) + m_trans[3];
	}

	/**
	 * The number of elements in the grid.
	 */
	size_t size() const {
		return (size_t) (m_cols * m_rows);
	}

	/**
	 * Returns true if the pixel is nodata.
	 */
	bool isNoData(int col, int row) {
		return get(col, row) == m_nodata;
	}

	bool isNoData(size_t idx) {
		return get(idx) == m_nodata;
	}

	/**
	 * Returns true if the pixel is nodata.
	 */
	bool isNoData(double x, double y) {
		return isNoData(toCol(x), toRow(y));
	}

	/**
	 * Returns true if the pixel exists and is not nodata.
	 */
	bool isValid(int c, int r) {
		return getOrNodata(c, r) != m_nodata;
	}

	/**
	 * Returns true if the pixel exists and is not nodata.
	 */
	bool isValid(double x, double y) {
		return getOrNodata(x, y) != m_nodata;
	}

	/**
	 * Gets the pixel value or nodata if the pixel doesn't exist.
	 */
	T getOrNodata(double x, double y) {
		if(!has(x, y)) {
			return m_nodata;
		} else {
			return get(x, y);
		}
	}

	/**
	 * Gets the pixel value or nodata if the pixel doesn't exist.
	 */
	T getOrNodata(int col, int row) {
		if(!has(col, row)) {
			return m_nodata;
		} else {
			return get(col, row);
		}
	}

	T *grid() {
		g_implerr("grid() Not implemented in Raster.");
	}
	
	bool hasGrid() const {
		return false;
	}

	/**
	 * Returns pixel value at the given coordinate.
	 */
	T &get(double x, double y) {
		return get(toCol(x), toRow(y));
	}

	/**
	 * Returns the pixel value at the give row/column.
	 */
	T &get(int col, int row) {
		loadBlock(col, row, false);
		return m_block[m_cache.toBlockIdx(col, row)];
	}

	T &get(size_t idx) {
		if(idx >= size())
			g_argerr("Index out of bounds.");
		return get(idx % m_cols, (int) idx / m_cols);
	}

	/**
	 * Return the element at the given index.
	 */
	T &operator[](size_t idx) {
		return get(idx);
	}

	/**
	 * Sets the pixel value at the given row/column.
	 */
	void set(int col, int row, T v) {
		if(!m_writable)
			g_runerr("This raster is not writable.");
		//g_trace("set " << col << ", " << row << ", " << v);
		loadBlock(col, row, true);
		m_block[m_cache.toBlockIdx(col, row)] = v;
	}

	void set(size_t idx, T v) {
		if(idx >= size())
			g_argerr("Index out of bounds.");
		set(idx % m_cols, (int) idx / m_cols, v);
	}

	/**
	 * Sets the pixel value at the given coordinate.
	 */
	void set(double x, double y, T v) {
		set(toCol(x), toRow(y), v);
	}

	bool isSquare() const {
		return cols() == rows();
	}

	/**
	 * Returns true if the col/row are represented in the dataset.
	 */
	bool has(int col, int row) const {
		return col >= 0 && col < m_cols && row >= 0 && row < m_rows;
	}

	/**
	 * Returns true if the x/y are represented in the dataset.
	 */
	bool has(double x, double y) const {
		return has(toCol(x), toRow(y));
	}

	bool has(size_t idx) const {
		return idx < (size_t) (m_cols * m_rows);
	}

	/**
	 * Flush the current block to the dataset.
	 */
	void flush() {
		if(m_writable)
			m_cache.flush();
		/*	
		if(m_writable && m_dirty) {
			if(m_band->WriteBlock(m_curcol, m_currow, m_block) != CE_None)
				g_runerr("Flush error.");
			m_ds->FlushCache();
			m_dirty = false;
		}
		*/
	}

	~Raster() {
		m_cache.close();
		if(m_ds) // Probably not necessary.
			GDALClose(m_ds);
	}

};

/**
 * The radius is given with cells as the unit, but
 * can be rational. When determining which cells to
 * include in the calculation, any cell which partially
 * falls in the radius will be included.
 */
template <class T>
void Grid<T>::voidFillIDW(double radius, int count, double exp) {

	if(radius <= 0.0)
		throw std::invalid_argument("Radius must be larger than 0.");

	if(count <= 0)
		throw std::invalid_argument("Count must be larger than 0.");

	if(exp <= 0.0)
		throw std::invalid_argument("Exponent must be larger than 0.");

	MemRaster<T> tmp(cols(), rows());
	tmp.nodata(nodata());
	tmp.fill(nodata());
	
	for(int r = 0; r < rows(); ++r) {
		for(int c = 0; c < cols(); ++c) {

			if(get(c, r) != nodata())
				continue;

			double rad = radius;
			bool found = false;

			do {

				double d = g_sq(rad);
				double a = 0.0;
				double b = 0.0;
				int cnt = 0;

				for(int r0 = g_max(0, r - rad); r0 < g_min(rows(), r + rad + 1); ++r0) {
					for(int c0 = g_max(0, c - rad); c0 < g_min(cols(), c + rad + 1); ++c0) {
						double d0 = g_sq((double) c0 - c) + g_sq((double) r0 - r);
						if(d0 <= d && get(c0, r0) != nodata()) {
							double dp = 1.0 / std::pow(d0, exp);
							a += dp * get(c0, r0);
							b += dp;
							++cnt;
						}
					}
				}

				if(cnt >= count) {
					tmp.set(c, r, a / b);
					found = true;
					break;
				}

				rad += 1.0;

			} while(rad < g_min(cols(), rows()));

			if(!found)
				std::cerr << "WARNING: Pixel not filled at " << c << "," << r << ". Consider larger radius or smaller count." << std::endl;
		}
	}

	writeBlock(tmp);
}


} // rater

} // geotools

#endif /* INCLUDE_RASTER_HPP_ */
