/*
 * Raster.hpp
 *
 *  Created on: Jan 21, 2016
 *      Author: rob
 */

#ifndef INCLUDE_RASTER_HPP_
#define INCLUDE_RASTER_HPP_

#include <queue>

#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <eigen3/Eigen/Core>

#include "geotools.h"

class Cell {
public:
	int col;
	int row;
	Cell(int col, int row) {
		this->col = col;
		this->row = row;
	}
};

template <class T>
class Grid {
protected:
	T m_min = 0;
	T m_max = 0;
	T m_mean = 0;
	T m_stddev = 0;
	bool m_stats = false;

public:
	
	virtual int rows() const =0;

	virtual int cols() const =0;
	
	virtual unsigned long size() const =0;

	virtual void fill(const T value) =0;

	/**
	 * Return a pointer to an in-memory grid of the raster data.
	 * Throw an exception if this is not possible.
	 */
	virtual T *grid() =0;
	
	/**
	 * Return a reference to the value held at
	 * the given index in the grid.
	 * Not const because the get operation might imply (e.g.)
	 * a buffering operation in the subclass.
	 */
	virtual T &get(unsigned long idx) =0;

	virtual T &get(int col, int row) =0;

	virtual void set(int col, int row, const T value) =0;

	virtual void set(unsigned long idx, const T value) =0;

	/**
	 * Return the element at the given index.
	 */
	virtual T &operator[](unsigned long idx) =0;

	virtual bool isSquare() const =0;

	virtual T nodata() const =0;

	virtual void nodata(T nodata) =0;

	/**
	 * Load the contents of a single block into the given Grid instance.
	 */
	virtual void readBlock(int col, int row, int cols, int rows, Grid<T> &block) =0;

	virtual void readBlock(int col, int row, Grid<T> &block) =0;

	/**
	 * Write a part of the given block to the raster.
	 */
	virtual void writeBlock(int col, int row, int cols, int rows, Grid<T> &block) =0;

	virtual void writeBlock(int col, int row, Grid<T> &block) =0;

	/**
	 * Write the full block to the raster.
	 */
	virtual void writeBlock(Grid<T> &block) =0;

	/**
	 * Read the full block from the raster.
	 */
	virtual void readBlock(Grid<T> &block) =0;

	void computeStats() {
		T total = 0;
		T max = get(0);
		T min = get(0);
		int count = 0;
		for(unsigned int i = 0; i < size(); ++i) {
			T v = get(i);
			if(v != nodata()) {
				total += v;
				if(v > max) max = v;
				if(v < min) min = v;
				++count;
			}
		}
		m_mean = total / count;
		m_min = min;
		m_max = max;
		T variance = 0;
		for(unsigned int i = 0; i < size(); ++i) {
			T v = get(i);
			if(v != nodata())
				variance += _sq(v - m_mean);
		}
		m_stddev = sqrt(variance / count);
		m_stats = true;
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

	void floodFill(int col, int row, T target, T fill, Grid<T> *other = nullptr, T otherFill = 0) {

		// If other is provided, fill pixels in it too, rather than just the current raster.

		T value = get(col, row);
		if(value != target)
			return;

		std::queue<Cell> q;
		q.push(Cell(col, row));

		while(q.size()) {
			Cell cel = q.front();
			q.pop();

			// Scan out to the left, starting with current.
			for(int c = cel.col; c >= 0; --c) {
				value = get(c, cel.row);
				if(value == target) {
					set(c, cel.row, fill);
					if(other)
						other->set(c, cel.row, otherFill);;
					if(cel.row > 0) {
						value = get(c, cel.row - 1);
						if(value == target)
							q.push(Cell(c, cel.row - 1));
					}
					if(cel.row < rows() - 1) {
						value = get(c, cel.row + 1);
						if(value == target) 
							q.push(Cell(c, cel.row + 1));
					}
				} else {
					break;
				}
			}

			// Scan out to the right.
			for(int c = cel.col + 1; c < cols(); ++c) {
				value = get(c, cel.row);
				if(value == target) {
					set(c, cel.row, fill);
					if(other)
						other->set(c, cel.row, otherFill);
					if(cel.row > 0) {
						value = get(c, cel.row - 1);
						if(value == target)
							q.push(Cell(c, cel.row - 1));
					}
					if(cel.row < rows() - 1) {
						value = get(c, cel.row + 1);
						if(value == target)
							q.push(Cell(c, cel.row + 1));
					}
				} else {
					break;
				}
			}
		}
	}

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
			throw "This instance has not been initialized.";
	}

public:
	MemRaster() {
		m_grid = nullptr;
		m_cols = -1;
		m_rows = -1;
		m_item_dealloc = nullptr;
		m_nodata = 0; // TODO: Choose a nodata based on type?
	}

	MemRaster(int cols, int rows) : MemRaster() {
		init(cols, rows);
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

	int rows() const {
		return m_rows;
	}

	int cols() const {
		return m_cols;
	}

	unsigned long size() const {
		return (unsigned long) m_rows * m_cols;
	}

	/**
	 * Initialize with the given number of cols and rows.
	 * (Re)allocates memory for the internal grid.
	 */
	void init(int cols, int rows) {
		if(cols <= 0 || rows <= 0)
			throw "Invalid row or column count.";
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
		std::fill_n(m_grid, size(), value);
	}

	/**
	 * Return a reference to the value held at
	 * the given index in the grid.
	 */
	T &get(unsigned long idx) {
		checkInit();
		if(idx >= size())
			throw "Index out of bounds.";
		return m_grid[idx];
	}

	T &get(int col, int row) {
		unsigned long idx = (unsigned long) (row * m_cols + col);	
		return get(idx);
	}

	void set(int col, int row, const T value) {
		unsigned long idx = (unsigned long) (row * m_cols + col);
		set(idx, value);
	}

	void set(unsigned long idx, const T value) {
		if(idx >= size())
			throw "Index out of bounds.";
		checkInit();
		m_grid[idx] = value;
	}

	/**
	 * Return the element at the given index.
	 */
	T &operator[](unsigned long idx) {
		checkInit();
		if(idx >= size())
			throw "Index out of bounds.";
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
	void readBlock(int col, int row, int cols, int rows, Grid<T> &block) {
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c)
				block.set(c, r, get(c + col, r + row));
		}
	}

	void readBlock(int col, int row, Grid<T> &block) {
		readBlock(col, row, _min(cols() - col, block.cols()), _min(rows() - row, block.rows()), block);
	}

	/**
	 * Write a part of the given block to the raster.
	 */
	void writeBlock(int col, int row, int cols, int rows, Grid<T> &block) {
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c)
				set(c + col, r + row, block.get(c, r));
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
	 * Read the full block from the raster.
	 */
	void readBlock(Grid<T> &block) {
		readBlock(0, 0, _min(block.cols(), cols()), _min(block.rows(), rows()), block);
	}


};

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
	T *m_block;				// Block storage
	GDALDataset *m_ds;		// GDAL dataset
	GDALRasterBand *m_band;	// GDAL band
	double m_trans[6];		// Raster transform
	bool m_inited = false;	// True if the instance is initialized.
	GDALDataType m_type;	// GDALDataType -- limits the possible template types.

	/**
	 * Loads the block that contains the given row and column.
	 */
	void loadBlock(int col, int row) {
		if(!has(col, row))
			throw "Row or column out of bounds.";
		int bcol = (int) (col / m_bw);
		int brow = (int) (row / m_bh);
		if(bcol >= m_bcols || bcol < 0 || brow >= m_brows || brow < 0)
			throw "Illegal block column or row.";
		if(bcol != m_curcol || brow != m_currow) {
			flush();
			if(m_band->ReadBlock(bcol, brow, m_block) != CE_None)
				throw "Failed to read block.";
			m_currow = brow;
			m_curcol = bcol;
		}
	}

	/**
	 * Flush the current block to the dataset.
	 */
	void flush() {
		if(m_writable && m_dirty) {
			if(m_band->WriteBlock(m_curcol, m_currow, m_block) != CE_None)
				throw "Flush error.";
			m_dirty = false;
		}
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
		m_writable(false), m_dirty(false) {
		m_ds = nullptr;
		m_band = nullptr;
		m_block = nullptr;
		m_type = getType((T) 0);
	}

	/**
	 * Create a new raster for writing with a template.
	 */
	template <class D>
	Raster(const std::string &filename, Raster<D> &tpl) : Raster() {
		std::string proj;
		tpl.projection(proj);
		init(filename, tpl.minx(), tpl.miny(), tpl.maxx(), tpl.maxy(), tpl.resolution(), (T) tpl.nodata(), proj);
	}

	Raster(const std::string &filename, Raster<T> &tpl) : Raster() {
		std::string proj;
		tpl.projection(proj);
		init(filename, tpl.minx(), tpl.miny(), tpl.maxx(), tpl.maxy(), tpl.resolution(), tpl.nodata(), proj);
	}

	/**
	 * Build a new raster with the given filename, bounds, resolution, nodata and projection.
	 */
	Raster(const std::string &filename, double minx, double miny, double maxx, double maxy,
			double resolution, double nodata, std::string &proj) : Raster() {
		init(filename, minx, miny, maxx, maxy, resolution, nodata, proj);
	}

	Raster(const std::string &filename, double minx, double miny, double maxx, double maxy,
			double resolution, double nodata, int crs) : Raster() {
		std::string proj = epsg2ProjText(crs);
		init(filename, minx, miny, maxx, maxy, resolution, nodata, proj);
	}

	/**
	 * Initializes the raster with the given filename, bounds, resolution, nodata and projection.
	 */
	void init(const std::string &filename, double minx, double miny, double maxx, double maxy,
			double resolution, double nodata, std::string proj) {
		if(minx >= maxx || miny >= maxy)
			throw "Invalid boundaries.";
		if(resolution <= 0.0)
			throw "Invalid resolution.";
		int width = (int) ((maxx - minx) / resolution) + 1;
		int height = (int) ((maxy - miny) / resolution) + 1;
		GDALAllRegister();
		m_ds = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(filename.c_str(),
				width, height, 1, m_type, NULL);
		if(m_ds == nullptr)
			throw "Failed to create file.";
		// TODO: Proper type.
		m_trans[0] = minx, m_trans[1] = resolution, m_trans[2] = 0.0,
				m_trans[3] = maxy, m_trans[4] = 0.0, m_trans[5] = -resolution;
		m_ds->SetGeoTransform(m_trans);
		if(!proj.empty())
			m_ds->SetProjection(proj.c_str());
		m_rows = m_ds->GetRasterYSize();
		m_cols = m_ds->GetRasterXSize();
		m_band = m_ds->GetRasterBand(1);
		if(m_band == NULL)
			throw "Failed to get band.";
		m_band->GetBlockSize(&m_bw, &m_bh);
		m_band->SetNoDataValue(nodata);
		m_nodata = m_band->GetNoDataValue();
		m_bcols = (m_cols + m_bw - 1) / m_bw;
		m_brows = (m_rows + m_bh - 1) / m_bh;
		m_block = (T *) malloc(sizeof(T) * m_bw * m_bh);
		m_writable = true;
		m_inited = true;
	}

	/**
	 * Open the given raster and load the given band. Set the writable argument to true
	 * to enable writing.
	 */
	Raster(std::string &filename, int band = 1, bool writable = false) : Raster() {
		init(filename, band, writable);
	}

	void init(std::string &filename, int band = 1, bool writable = false) {
		GDALAllRegister();
		m_ds = (GDALDataset *) GDALOpen(filename.c_str(), writable ? GA_Update : GA_ReadOnly);
		if(m_ds == NULL)
			throw "Failed to open raster.";
		m_bandn = band;
		m_ds->GetGeoTransform(m_trans);
		m_band = m_ds->GetRasterBand(band);
		if(m_band == nullptr)
			throw "Failed to get band.";
		m_band->GetBlockSize(&m_bw, &m_bh);
		m_rows = m_ds->GetRasterYSize();
		m_cols = m_ds->GetRasterXSize();
		m_nodata = m_band->GetNoDataValue();
		m_bcols = (m_cols + m_bw - 1) / m_bw;
		m_brows = (m_rows + m_bh - 1) / m_bh;
		m_block = (T *) malloc(sizeof(T) * m_bw * m_bh);
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
		if(m_band->RasterIO(GF_Read, col, row, cols, rows, block.grid(), cols, rows, m_type, 0, 0) != CE_None)
			throw "Error reading block.";
	}

	void readBlock(int col, int row, Grid<T> &block) {
		readBlock(col, row, _min(cols() - col, block.cols()), _min(rows() - row, block.rows()), block);
	}

	/**
	 * Write a part of the given block to the raster.
	 */
	void writeBlock(int col, int row, int cols, int rows, Grid<T> &block) {
		if(m_band->RasterIO(GF_Write, col, row, cols, rows, block.grid(), cols, rows, m_type, 0, 0) != CE_None)
			throw "Error writing block.";
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
	 * Read the full block from the raster.
	 */
	void readBlock(Grid<T> &block) {
		readBlock(0, 0, _min(block.cols(), cols()), _min(block.rows(), rows()), block);
	}
	
	/**
	 * Return the resolution. If the x and y resolution are different,
	 * use resolutionX and resolutionY.
	 */
	double resolution() const {
		return m_trans[1];
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
		throw "Cannot set nodata on Raster.";
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

	/**
	 * The number of elements in the grid.
	 */
	unsigned long size() const {
		return (unsigned long) (m_cols * m_rows);
	}

	/**
	 * Returns true if the pixel is nodata.
	 */
	bool isNoData(int col, int row) {
		return get(col, row) == m_nodata;
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
		throw "grid() Not implemented in Raster.";
	}
	
	/**
	 * Returns pixel value at the given coordinate.
	 */
	T &get(double x, double y) const {
		return get(toCol(x), toRow(y));
	}

	/**
	 * Returns the pixel value at the give row/column.
	 */
	T &get(int col, int row) {
		loadBlock(col, row);
		unsigned long idx = (unsigned long) (row % m_bh) * m_bw + (col % m_bw);
		return m_block[idx];
	}

	T &get(unsigned long idx) {
		int col = idx % m_bw;
		int row = (int) (idx / m_bw);
		return get(col, row);
	}

	/**
	 * Return the element at the given index.
	 */
	T &operator[](unsigned long idx) {
		int col = idx % m_bw;
		int row = (int) (idx / m_bw);
		return get(col, row);
	}

	/**
	 * Sets the pixel value at the given row/column.
	 */
	void set(int col, int row, T v) {
		if(!m_writable) return;
		loadBlock(col, row);
		unsigned long idx = (unsigned long) (row % m_bh) * m_bw + (col % m_bw);
		m_block[idx] = v;
		m_dirty = true;
	}

	void set(unsigned long idx, T v) {
		if(idx >= size())
			throw "Index out of bounds.";
		m_block[idx] = v;
		m_dirty = true;
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

	~Raster() {
		flush();
	}

};

#endif /* INCLUDE_RASTER_HPP_ */
