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
#include <memory>

#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <eigen3/Eigen/Core>

#include "geotools.h"

namespace raster {

namespace util {

	/**
	 * These methods return the GDAL data type given the template
	 * type of the class.
	 */
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

}

/**
 * Convenience class for storing row/column coordinates.
 */
class Cell {
public:
	int col;
	int row;
	/**
	 * Construct an instance using the given row and column.
	 */
	Cell(int col, int row) {
		this->col = col;
		this->row = row;
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

/**
 * A simple FillOperator that simply matches the given value
 * to the given target.
 */
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

public:
	
	virtual int rows() const =0;
	virtual int cols() const =0;
	
	/**
	 * The number of pixels in the raster.
	 */
	virtual unsigned long size() const =0;

	/**
	 * Fill the entire raster with the given value.
	 */
	virtual void fill(const T value) =0;

	/**
	 * Return a pointer to an in-memory grid of the raster data.
	 * Throw an exception if this is not possible.
	 */
	virtual T *grid() =0;
	
	/**
	 * Return a reference to the value held at the given index in the grid.
	 * Not const because the get operation might imply (e.g.) a buffering operation in the subclass.
	 */
	virtual T &get(unsigned long idx) =0;
	virtual T &get(int col, int row) =0;
	virtual T &operator[](unsigned long idx) =0;

	/**
	 * Set the value at the given position in the raster. 
	 * Not const because the get operation might imply (e.g.) a buffering operation in the subclass.
	 */
	virtual void set(int col, int row, const T value) =0;
	virtual void set(unsigned long idx, const T value) =0;

	virtual bool has(int col, int row) const =0;
	virtual bool has(unsigned long idx) const =0;

	/**
	 * Returns true if the value at the given position is a nodata value.
 	 * Not const because the get operation might imply (e.g.) a buffering operation in the subclass.
	 */
	virtual bool isNoData(int col, int row) =0;
	virtual bool isNoData(unsigned long idx) =0;

	/**
	 * Returns true if a nodata value has been explicitly set.
	 */
	virtual bool hasNodata() const =0;

	/**
	 * Returns true if the raster is a square.
	 */
	virtual bool isSquare() const =0;

	virtual T nodata() const =0;
	virtual void nodata(T nodata) =0;

	/**
	 * Load the contents of a single block into the given Grid instance.
	 * col and row parameters give the location of the start of the read in the source raster.
	 * cols and rows give the width and height of the read block.
	 */
	virtual void readBlock(int col, int row, int cols, int rows, Grid<T> &block) =0;
	virtual void readBlock(int col, int row, Grid<T> &block) =0;
	virtual void readBlock(Grid<T> &block) =0;

	/**
	 * Write a part of the given block to the raster.
	 * col and row parameters give the location of the start of the write in the destination raster.
	 * cols and rows give the width and height of the written block.
	 */
	virtual void writeBlock(int col, int row, int cols, int rows, Grid<T> &block) =0;
	virtual void writeBlock(int col, int row, Grid<T> &block) =0;
	virtual void writeBlock(Grid<T> &block) =0;

	/**
	 * Return a "slice" of the raster.
	 * col and row are the location of the start of reading in the source raster.
	 * cols and rows are the width and height of the read region.
	 * A new Grid is constructed and returned with the same dimensions given.
	 */
	virtual std::unique_ptr<Grid<T> > slice(int col, int row, int cols, int rows) =0;

	/**
	 * Compute statistics for this grid. 
	 * Statistics are, min, max, mean, variance and std. deviation.
	 * A nodata value must be set if the dataset has any such values, otherwise
	 * these statistics are bogus. Therefore an exception is thrown is it is not set.
	 */
	void computeStats() {
		if(!hasNodata())
			_runerr("No nodata value set. Statistics cannot be computed.");
		
		// Locate the first non-nodata value.
		unsigned int i;
		for(i = 0; i < size(); ++i) {
			if(!isNoData(i)) {
				m_min = m_max = get(i);
				break;
			}
		}

		// Welford's method for variance.
		// i has the index of the first non-nodata element.
		m_sum = 0;
		m_count = 0;
		T m = 0;
		T s = 0;
		int k = 1;
		for(; i < size(); ++i) {
			T v = get(i);
			if(v != nodata()) {
				T oldm = m;
				m = m + (v - m) / k;
				s = s + (v - m) * (v - oldm);
				m_sum += v;
				m_min = _min(m_min, v);
				m_max = _max(m_max, v);
				++m_count;
				++k;
			}
		}
		m_mean = m_sum / m_count;
		m_variance = s / m_count;
		m_stddev = std::sqrt(m_variance);
		m_stats = true;
		_trace("Count: " << m_count << "; Sum: " << m_sum << "; Min: " << m_min << "; Max: " 
			<< m_max << "; Mean: " << m_mean 
			<< "; Variance: " << m_variance << "; Std Dev: " << m_stddev << std::endl);
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
	 * Performs a flood fill operation beginning at the given column and row.
	 * The FillOperator determines whether a target is valid. 
	 * The fill is actually applied to the other Grid, using the current grid as the target for
	 * matching purposes.
	 */
	template <class U>
	std::vector<int> floodFill(int col, int row, FillOperator<T> &op, Grid<U> &other, U otherFill) {

		int minc = cols() + 1;
		int minr = rows() + 1;
		int maxc = -1;
		int maxr = -1;
		int area = 0;

		if(other.get(col, row) == otherFill || !op.fill(get(col, row)))
			return {minc, minr, maxc, maxr, area};

		std::queue<Cell> q;
		q.push(Cell(col, row));

		while(q.size()) {
			Cell cel = q.front();
			q.pop();

			other.set(cel.col, cel.row, otherFill);

			minc = _min(cel.col, minc);
			maxc = _max(cel.col, maxc);
			minr = _min(cel.row, minr);
			maxr = _max(cel.row, maxr);
			++area;

			// Scan out to the left, starting with current.
			for(int c = cel.col - 1; c >= 0; --c) {
				if(other.get(c, cel.row) != otherFill && op.fill(get(c, cel.row))) {

					minc = _min(c, minc);
					++area;

					other.set(c, cel.row, otherFill);

					if(cel.row > 0) {
						if(other.get(c, cel.row - 1) != otherFill && op.fill(get(c, cel.row - 1)))
							q.push(Cell(c, cel.row - 1));
					}
					if(cel.row < rows() - 1) {
						if(other.get(c, cel.row + 1) != otherFill && op.fill(get(c, cel.row + 1)))
							q.push(Cell(c, cel.row + 1));
					}
				} else {
					break;
				}
			}

			// Scan out to the right.
			for(int c = cel.col + 1; c < cols(); ++c) {
				if(other.get(c, cel.row) != otherFill && op.fill(get(c, cel.row))) {

					maxc = _max(c, maxc);
					++area;

					other.set(c, cel.row, otherFill);

					if(cel.row > 0) {
						if(other.get(c, cel.row - 1) != otherFill && op.fill(get(c, cel.row - 1)))
							q.push(Cell(c, cel.row - 1));
					}
					if(cel.row < rows() - 1) {
						if(other.get(c, cel.row + 1) != otherFill && op.fill(get(c, cel.row + 1)))
							q.push(Cell(c, cel.row + 1));
					}
				} else {
					break;
				}
			}
		}
		return {minc, minr, maxc, maxr, area};
	}

	/**
	 * Fill the present grid with the fill using the given target.
	 */
	std::vector<int> floodFill(int col, int row, T target, T fill);

	/**
	 * Fill the present grid with the fill using the given operator.
	 */	
	std::vector<int> floodFill(int col, int row, FillOperator<T> &op, T fill);
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
	bool m_nodataSet;

	void (*m_item_dealloc)(T);

	/**
	 * Checks if the grid has been initialized. Throws exception otherwise.
	 */
	void checkInit() const {
		if(m_grid == nullptr)
			_runerr("This instance has not been initialized.");
	}

public:
	MemRaster() :
		m_cols(-1), m_rows(-1),
		m_grid(nullptr),
		m_nodata(0), m_nodataSet(false),
		m_item_dealloc(nullptr) {
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
			for(unsigned long i = 0; i < (unsigned long) m_cols * m_rows; ++i)
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

	/**
	 * Cast a MemRaster to some other type. Might cause bad things.
	 */
	template <class U>
	operator MemRaster<U>() {
		MemRaster<U> g(cols(), rows());
		for(unsigned long i = 0; i < size(); ++i)
			g.set(i, (U) get(i));
		return g;
	}

	int rows() const {
		return m_rows;
	}

	int cols() const {
		return m_cols;
	}

	unsigned long size() const {
		return (unsigned long) m_rows * m_cols;
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
			_argerr("Invalid row or column count.");
		if(cols != m_cols || rows != m_rows) {
			m_cols = cols;
			m_rows = rows;
			if(m_grid != nullptr)
				free(m_grid);
			m_grid = (T *) malloc(sizeof(T) * cols * rows);
		}
	}

	void fill(const T value) {
		checkInit();
		std::fill_n(m_grid, size(), value);
	}

	std::unique_ptr<Grid<T> > slice(int col, int row, int cols, int rows) {
		Grid<T> *g = new MemRaster<T>(cols, rows);
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c)
				g->set(c, r, get(c + col, r + row));
		}
		return std::unique_ptr<Grid<T> >(g);
	}

	T &get(unsigned long idx) {
		checkInit();
		if(idx >= size())
			_argerr("Index out of bounds: size: " << size() << "; idx: " << idx);
		return m_grid[idx];
	}

	T &get(int col, int row) {
		unsigned long idx = (unsigned long) (row * m_cols + col);	
		return get(idx);
	}

	T &operator[](unsigned long idx) {
		return get(idx);
	}

	bool isNoData(int col, int row) {
		return get(col, row) == m_nodata;
	}

	bool isNoData(unsigned long idx) {
		return get(idx) == m_nodata;
	}

	void set(int col, int row, const T value) {
		unsigned long idx = (unsigned long) (row * m_cols + col);
		set(idx, value);
	}

	void set(unsigned long idx, const T value) {
		if(idx >= size())
			_argerr("Index out of bounds: size: " << size() << "; idx: " << idx);
		checkInit();
		m_grid[idx] = value;
	}

	bool has(int col, int row) const {
		return col >= 0 && col < m_cols && row >= 0 && row < m_rows;
	}

	bool has(unsigned long idx) const {
		return idx < (unsigned long) (m_cols * m_rows);
	}

	bool isSquare() const {
		return cols() == rows();
	}

	/**
	 * Configures an Eigen::Matrix with the contents of this grid.
	 */
	void toMatrix(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx) {
		for(int r = 1; r < rows(); ++r) {
			for(int c = 0; c < cols(); ++c)
				mtx(r, c) = get(c, r);
		}
	}

	/**
	 * Configures this grid with the contents of the Eigen::Matrix.
	 */
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
		m_nodataSet = true;
	}

	bool hasNodata() const {
		return m_nodataSet;
	}

	void readBlock(int col, int row, int cols, int rows, Grid<T> &block) {
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c)
				block.set(c, r, get(c + col, r + row));
		}
	}

	void readBlock(int col, int row, Grid<T> &block) {
		readBlock(col, row, _min(cols() - col, block.cols()), _min(rows() - row, block.rows()), block);
	}

	void readBlock(Grid<T> &block) {
		readBlock(0, 0, _min(block.cols(), cols()), _min(block.rows(), rows()), block);
	}

	void writeBlock(int col, int row, int cols, int rows, Grid<T> &block) {
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c)
				set(c + col, r + row, block.get(c, r));
		}
	}

	void writeBlock(int col, int row, Grid<T> &block) {
		writeBlock(col, row, _min(cols() - col, block.cols()), _min(rows() - row, block.rows()), block);
	}

	void writeBlock(Grid<T> &block) {
		writeBlock(0, 0, _min(block.cols(), cols()), _min(block.rows(), rows()), block);
	}


};

/**
 * This is a helper class that permits the user to iterate over the blocks in a
 * Raster.
 */
template <class T>
class Block {
private:
	int m_bc, m_br;
	int m_bw, m_bh;
	int m_col, m_row;
	int m_tc, m_tr;

public:
	Block(int blockCols, int blockRows, int blockWidth, int blockHeight, int totalCols, int totalRows) {
		m_bc = blockCols;
		m_br = blockRows;
		m_bw = blockWidth;
		m_bh = blockHeight;
		m_tc = totalCols;
		m_tr = totalRows;
		m_col = -1;
		m_row = 0;
	}

	bool next() {
		if(++m_col == m_bc) {
			m_col = 0;
			++m_row;
		}
		return m_col < m_bc && m_row < m_br;
	}
	void reset() {
		m_col = m_row = 0;
	}
	int rows() {
		return endRow() - startRow();
	}
	int cols() {
		return endCol() - startCol();
	}
	int startCol() {
		return m_col * m_bw;
	}
	int endCol() {
		return _min(startCol() + m_bw, m_tc);
	}
	int startRow() {
		return m_row * m_bh;
	}
	int endRow() {
		return _min(startRow() + m_bh, m_tr);
	}
	~Block() {}
};

/**
 * This class represents a file-based raster.
 */
template <class T>
class Raster : public Grid<T> {
private:
	int m_cols, m_rows;		// Raster cols/rows
	int m_bcols, m_brows;	// Block cols/rows -- not the number of cols/rows in a block
	int m_curcol, m_currow;	// The current block column
	int m_bandn;			// The band number
	int m_bw, m_bh;			// Block width/height in pixels
	bool m_writable;		// True if the raster is writable
	bool m_dirty;			// True if there is a modification that should be flushed.
	T m_nodata;				// Nodata value.
	bool m_nodataSet;       // True if nodata is explicitly set.
	T *m_block;				// Block storage
	GDALDataset *m_ds;		// GDAL dataset
	GDALRasterBand *m_band;	// GDAL band
	double m_trans[6];		// Raster transform
	bool m_inited;	        // True if the instance is initialized.
	GDALDataType m_type;	// GDALDataType -- limits the possible template types.

	/**
	 * Loads the block that contains the given row and column.
	 */
	void loadBlock(int col, int row) {
		if(!m_inited)
			_runerr("Not inited before attempted read.");
		if(!has(col, row))
			_argerr("Row or column out of bounds.");
		int bcol = (int) (col / m_bw);
		int brow = (int) (row / m_bh);
		if(bcol >= m_bcols || bcol < 0 || brow >= m_brows || brow < 0)
			_argerr("Illegal block column or row.");
		if(bcol != m_curcol || brow != m_currow) {
			flush();
			if(m_band->ReadBlock(bcol, brow, m_block) != CE_None)
				_runerr("Failed to read block.");
			m_currow = brow;
			m_curcol = bcol;
		}
	}

	/**
	 * Flush the current block to the dataset.
	 */
	void flush() {
		if(m_writable && m_dirty) {
			std::cerr << "write " << m_curcol << " " << m_currow << " " << m_block << std::endl;
			if(m_band->WriteBlock(m_curcol, m_currow, m_block) != CE_None)
				_runerr("Flush error.");
			m_dirty = false;
		}
	}

public:

	/**
	 * Basic constructor.
	 */
	Raster() :
		m_cols(-1), m_rows(-1),
		m_bcols(-1), m_brows(-1),
		m_curcol(-1), m_currow(-1),
		m_bandn(1),
		m_bw(-1), m_bh(-1),
		m_writable(false), m_dirty(false),
		m_nodata(0), m_nodataSet(false),
		m_block(nullptr),
		m_ds(nullptr), m_band(nullptr),
		m_inited(false),
		m_type(util::getType((T) 0)) {
	}

	/**
	 * Create a new raster for writing with a template.
	 */
	template <class D>
	Raster(const std::string &filename, Raster<D> &tpl) : Raster() {
		std::string proj;
		tpl.projection(proj);
		init(filename, tpl.minx(), tpl.miny(), tpl.maxx(), tpl.maxy(), tpl.resolutionX(), 
			tpl.resolutionY(), (T) tpl.nodata(), proj);
	}

	/**
	 * Build a new raster with the given filename, bounds, resolution, nodata and projection.
	 */
	Raster(const std::string &filename, double minx, double miny, double maxx, double maxy,
			double resolutionX, double resolutionY, double nodata, std::string &proj) : Raster() {
		init(filename, minx, miny, maxx, maxy, resolutionX, resolutionY, nodata, proj);
	}

	/**
	 * Build a new raster with the given filename, bounds, resolution, nodata and SRID.
	 */
	Raster(const std::string &filename, double minx, double miny, double maxx, double maxy,
			double resolutionX, double resolutionY, double nodata, int crs) : Raster() {
		std::string proj = epsg2ProjText(crs);
		init(filename, minx, miny, maxx, maxy, resolutionX, resolutionY, nodata, proj);
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
	void init(const std::string &filename, Raster<D> &tpl) {
		std::string proj;
		tpl.projection(proj);
		init(filename, tpl.minx(), tpl.miny(), tpl.maxx(), tpl.maxy(), tpl.resolutionX(), 
			tpl.resolutionY(), tpl.nodata(), proj);
	}

	/**
	 * Initializes the raster with the given filename, bounds, resolution, nodata and projection.
	 */
	void init(const std::string &filename, double minx, double miny, double maxx, double maxy,
			double resolutionX, double resolutionY, double nodata, std::string proj) {
		
		if(resolutionX == 0 || resolutionY == 0)
			_argerr("Resolution must be larger or smaller than 0.");
		if(maxx < minx)
			_argerr("Minimum x must be smaller than or equal to maximum x.");
		if(maxy < miny)
			_argerr("Minimum y must be smaller than or equal to maximum y.");
		
		// Compute columns/rows
		int width = (int) ((maxx - minx) / resolutionX) * (resolutionX < 0 ? -1 : 1);
		int height = (int) ((maxy - miny) / resolutionY) * (resolutionY < 0 ? -1 : 1);

		// Create GDAL dataset.
		GDALAllRegister();
		m_ds = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(filename.c_str(),
				width, height, 1, m_type, NULL);
		if(m_ds == nullptr)
			_runerr("Failed to create file.");

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
		m_band = m_ds->GetRasterBand(1);
		if(m_band == NULL)
			_runerr("Failed to get band.");
		m_band->GetBlockSize(&m_bw, &m_bh);
		m_band->SetNoDataValue(nodata);
		m_nodata = m_band->GetNoDataValue();
		m_nodataSet = true;
		m_bcols = m_cols / m_bw;
		m_brows = m_rows / m_bh;
		m_block = (T *) malloc(sizeof(T) * m_bw * m_bh);
		if(!m_block)
			_runerr("Failed to allocate memory for raster block.");
		m_writable = true;
		m_inited = true;
	}

	/**
	 * Initializes a Raster from the existing file.
	 */
	void init(const std::string &filename, int band = 1, bool writable = false) {
		// Attempt to open the dataset.
		GDALAllRegister();
		m_ds = (GDALDataset *) GDALOpen(filename.c_str(), writable ? GA_Update : GA_ReadOnly);
		if(m_ds == NULL)
			_runerr("Failed to open raster.");

		// Save some raster
		m_bandn = band;
		m_ds->GetGeoTransform(m_trans);
		m_band = m_ds->GetRasterBand(band);
		if(m_band == nullptr)
			_runerr("Failed to get band.");
		m_band->GetBlockSize(&m_bw, &m_bh);
		m_rows = m_ds->GetRasterYSize();
		m_cols = m_ds->GetRasterXSize();
		m_nodata = m_band->GetNoDataValue();
		m_nodataSet = true;
		m_bcols = m_cols / m_bw;
		m_brows = m_rows / m_bh;
		m_block = (T *) malloc(sizeof(T) * m_bw * m_bh);
		if(!m_block)
			_runerr("Failed to allocate memory for raster block.");
		m_writable = writable;
		m_inited = true;
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
		Block<T> blk = block();
		MemRaster<T> grd(blk.cols(), blk.rows());
		grd.fill(value);
		while(blk.next())
			writeBlock(blk.startCol(), blk.startRow(), grd);
	}

	/**
	 * Return a "Block" which just stores indices of pixels that
	 * comprise a block. Calling next on the block adjusts the indices
	 * to correspond to the next block.
	 */
	Block<T> block() {
		return Block<T>(m_bcols, m_brows, m_bw, m_bh, m_cols, m_rows);
		// https://en.wikipedia.org/wiki/Return_value_optimization
	}

	/**
	 * Return the number of blocks in this raster.
	 */
	unsigned int numBlocks() {
		return (unsigned int) m_bcols * m_brows;
	}

	/**
	 * Load the contents of a single block into the given Grid instance.
	 */
	void readBlock(int col, int row, int cols, int rows, Grid<T> &block) {
		// TODO: Use readblock/writeblock
		if(col % m_bw == 0 && row % m_bh == 0 && block.cols() == m_bw && block.rows() == m_bh) {
			if(m_band->ReadBlock(col / m_bw, row / m_bh, block.grid()) != CE_None)
				_runerr("Error reading block (1).");
		} else {
			if(m_band->RasterIO(GF_Read, col, row, cols, rows, block.grid(), cols, rows, m_type, 0, 0) != CE_None)
				_runerr("Error reading block (2).");
		}	
	}

	void readBlock(int col, int row, Grid<T> &block) {
		readBlock(col, row, _min(cols() - col, block.cols()), _min(rows() - row, block.rows()), block);
	}

	void readBlock(Grid<T> &block) {
		readBlock(0, 0, _min(block.cols(), cols()), _min(block.rows(), rows()), block);
	}

	/**
	 * Write a part of the given block to the raster.
	 */
	void writeBlock(int col, int row, int cols, int rows, Grid<T> &block) {
		if(col % m_bw == 0 && row % m_bh == 0 && block.cols() == m_bw && block.rows() == m_bh) {
			if(m_band->WriteBlock(col / m_bw, row / m_bh, block.grid()) != CE_None)
				_runerr("Error writing block (1).");
		} else {
			if(m_band->RasterIO(GF_Write, col, row, cols, rows, block.grid(), cols, rows, m_type, 0, 0) != CE_None)
				_runerr("Error writing block (2).");
		}
	}

	void writeBlock(int col, int row, Grid<T> &block) {
		writeBlock(col, row, _min(cols() - col, block.cols()), _min(rows() - row, block.rows()), block);
	}

	/**
	 * Write the full block to the raster.
	 */
	void writeBlock(Grid<T> &block) {
		writeBlock(0, 0, _min(block.cols(), cols()), _min(block.rows(), rows()), block);
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

	std::unique_ptr<Grid<T> > slice(int col, int row, int cols, int rows) {
		Grid<T> *g = new MemRaster<T>(cols, rows);
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c)
				g->set(c, r, get(c + col, r + row));
		}
		return std::unique_ptr<Grid<T> >(g);
	}

	/**
	 * Return the GDAL datatype of the raster.
	 */
	GDALDataType type() const {
		return m_band->GetRasterDataType();
	}

	/**
	 * Return the leftmost X coordinate of the raster. For negative resolution, 
	 * this could be larger than the rightmost.
	 */
	double lx() {
		return m_trans[0];
	}

	/**
	 * Return the rightmost X coordinate of the raster. For negative resolution, 
	 * this could be smaller than the leftmost.
	 */
	double rx() {
		return m_trans[0] + m_cols * m_trans[1];
	}

	/**
	 * Return the top X coordinate of the raster. For negative resolution, 
	 * this could be larger than the bottom.
	 */
	double ty() {
		return m_trans[3];
	}

	/**
	 * Return the bottom Y coordinate of the raster. For negative resolution, 
	 * this could be smaller than the top.
	 */
	double by() {
		return m_trans[3] + m_rows * m_trans[5];
	}

	/**
	 * The minimum bounding x.
	 */
	double minx() const {
		return m_trans[1] > 0 ? m_trans[0] : m_trans[0] + m_cols * m_trans[1];
	}

	/**
	 * The maximum bounding x.
	 */
	double maxx() const {
		return m_trans[1] < 0 ? m_trans[0] : m_trans[0] + m_cols * m_trans[1];
	}

	/**
	 * The minimum bounding y.
	 */
	double miny() const {
		return m_trans[5] > 0 ? m_trans[3] : m_trans[3] + m_cols * m_trans[5];
	}

	/**
	 * The maximum bounding y.
	 */
	double maxy() const {
		return m_trans[5] > 0 ? m_trans[3] : m_trans[3] + m_cols * m_trans[5];
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

	T nodata() const {
		return m_nodata;
	}

	void nodata(T nodata) {
		try {
			m_band->SetNoDataValue(nodata);
			m_nodataSet = true;
		} catch(...) {
			_runerr("Cannot set nodata on Raster.");
		}
	}

	bool hasNodata() const {
		return m_nodataSet;
	}

	/*
	 * Returns the row offset in the block for a given y.
	 */
	int toBlockRow(double y) const {
		return toRow(y) % m_brows;
	}

	/**
	 * Returns the row offset in the block for a given x.
	 */
	int toBlockCol(double x) const {
		return toCol(x) % m_bcols;
	}

	/**
	 * Returns the width of the block (number of cells).
	 */
	int blockWidth() const {
		return m_bw;
	}

	/**
	 * Returns the height of the block (number of cells).
	 */
	int blockHeight() const {
		return m_bh;
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

	unsigned long size() const {
		return (unsigned long) (m_cols * m_rows);
	}

	/**
	 * Returns true if the pixel is nodata.
	 */
	bool isNoData(int col, int row) {
		return get(col, row) == m_nodata;
	}

	bool isNoData(unsigned long idx) {
		return get(idx) == m_nodata;
	}

	bool isNoData(double x, double y) {
		return isNoData(toCol(x), toRow(y));
	}

	/**
	 * Returns true if the pixel exists and is not nodata.
	 */
	bool isValid(int c, int r) {
		return getOrNodata(c, r) != m_nodata;
	}

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

	T getOrNodata(int col, int row) {
		if(!has(col, row)) {
			return m_nodata;
		} else {
			return get(col, row);
		}
	}

	T *grid() {
		_implerr("grid() Not implemented in Raster.");
	}
	
	T &get(double x, double y) const {
		return get(toCol(x), toRow(y));
	}

	T &get(int col, int row) {
		loadBlock(col, row);
		unsigned long idx = (unsigned long) (row % m_bh) * m_bw + (col % m_bw);
		return m_block[idx];
	}

	T &get(unsigned long idx) {
		if(idx >= size())
			_argerr("Index out of bounds: size: " << size() << "; idx: " << idx);
		return get(idx % m_cols, (int) idx / m_rows);
	}

	T &operator[](unsigned long idx) {
		return get(idx);
	}

	void set(int col, int row, T v) {
		if(!m_writable) return;
		loadBlock(col, row);
		unsigned long idx = (unsigned long) (row % m_bh) * m_bw + (col % m_bw);
		m_block[idx] = v;
		m_dirty = true;
	}

	void set(unsigned long idx, T v) {
		if(idx >= size())
			_argerr("Index out of bounds: size: " << size() << "; idx: " << idx);
		set(idx % m_cols, (int) idx / m_rows, v);
	}

	void set(double x, double y, T v) {
		set(toCol(x), toRow(y), v);
	}

	bool isSquare() const {
		return cols() == rows();
	}

	bool has(int col, int row) const {
		return col >= 0 && col < m_cols && row >= 0 && row < m_rows;
	}

	bool has(double x, double y) const {
		return has(toCol(x), toRow(y));
	}

	bool has(unsigned long idx) const {
		return idx < (unsigned long) (m_cols * m_rows);
	}

	~Raster() {
		flush();
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
				double d = _sq(rad);
				double a = 0.0;
				double b = 0.0;
				int cnt = 0;
				for(int r0 = _max(0, r - rad); r0 < _min(rows(), r + rad + 1); ++r0) {
					for(int c0 = _max(0, c - rad); c0 < _min(cols(), c + rad + 1); ++c0) {
						double d0 = _sq((double) c0 - c) + _sq((double) r0 - r);
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

			} while(rad < _min(cols(), rows()));

			if(!found)
				_warn("Pixel not filled at " << c << "," << r << ". Consider larger radius or smaller count.");
		}
	}

	writeBlock(tmp);
}

/**
 * Fill the present grid with the fill using the given target.
 */
template <class T>
std::vector<int> Grid<T>::floodFill(int col, int row, T target, T fill) {
	TargetOperator<T> op(target);
	MemRaster<T> tmp(*this);
	return floodFill(col, row, op, tmp, fill);
}

/**
 * Fill the present grid with the fill using the given operator.
 */	
template <class T>
std::vector<int> Grid<T>::floodFill(int col, int row, FillOperator<T> &op, T fill) {
	MemRaster<T> tmp(*this);
	return floodFill(col, row, op, tmp, fill);

}


} // raster

#endif /* INCLUDE_RASTER_HPP_ */
