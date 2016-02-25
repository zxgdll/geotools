/*
 * Raster.hpp
 *
 *  Created on: Jan 21, 2016
 *      Author: rob
 */

#ifndef INCLUDE_RASTER_HPP_
#define INCLUDE_RASTER_HPP_

#include <gdal_priv.h>
#include <ogr_spatialref.h>

template <class T>
class Block {
private:
	int m_bc, m_br;
	int m_bw, m_bh;
	int m_col, m_row;
	int m_tc, m_tr;
	int min(int a, int b) {
		return a < b ? a : b;
	}
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
	int startCol() {
		return m_col * m_bw;
	}
	int endCol() {
		return min(startCol() + m_bw, m_tc);
	}
	int startRow() {
		return m_row * m_bh;
	}
	int endRow() {
		return min(startRow() + m_bh, m_tr);
	}
	~Block() {}
};

template <class T>
class Raster {
private:
	int m_cols, m_rows;			// Raster cols/rows
	int m_bcols, m_brows;		// Block cols/rows -- not the number of cols/rows in a block
	int m_curcol, m_currow;		// The current block column
	int m_bandn;				// The band number
	int m_bw, m_bh;				// Block width/height in pixels
	bool m_writable;			// True if the raster is writable
	bool m_dirty;
	T m_nodata;					// Nodata value.
	T *m_block;					// Block storage
	GDALDataset *m_ds;			// GDAL dataset
	GDALRasterBand *m_band;		// GDAL band
	double m_trans[6];			// Raster transform
	bool m_inited = false;

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
	}

	/**
	 * Build a new raster with the given filename, bounds, resolution, nodata and projection.
	 */
	Raster(std::string &filename, double minx, double miny, double maxx, double maxy,
			double resolution, double nodata, const char *proj = NULL) : Raster() {
			init(filename, minx, miny, maxx, maxy, resolution, nodata, proj);
	}

	void init(std::string &filename, double minx, double miny, double maxx, double maxy,
			double resolution, double nodata, const char *proj = NULL) {
		GDALAllRegister();
		int width = (int) ((maxx - minx) / resolution) + 1;
		int height = (int) ((maxy - miny) / resolution) + 1;
		m_ds = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(filename.c_str(),
				width, height, 1, GDT_Float32, NULL);
		if(m_ds == nullptr)
			throw "Failed to create file.";
		// TODO: Proper type.
		m_trans[0] = minx, m_trans[1] = resolution, m_trans[2] = 0.0,
				m_trans[3] = maxy, m_trans[4] = 0.0, m_trans[5] = -resolution;
		m_ds->SetGeoTransform(m_trans);
		if(proj != NULL)
			m_ds->SetProjection(proj);
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

	bool inited() {
		return m_inited;
	}
	/**
	 * Return a "Block" which just stores indices of pixels that
	 * comprise a block. Calling next on the block adjusts the indices
	 * to correspond to the next block.
	 */
	std::unique_ptr<Block<T> > block() {
		return std::unique_ptr<Block<T> >(new Block<T>(m_bcols, m_brows, m_bw, m_bh, m_cols, m_rows));
	}

	/**
	 * Return the number of blocks in this raster.
	 */
	int numBlocks() {
		return m_bcols * m_brows;
	}

	/**
	 * Load the contents of a single block into the given array.
	 * Caller is responsible for correctly initializing the array,
	 * and freeing it.
	 */
	void loadBlock(int col, int row, int cols, int rows, T *block) {
		// TODO: Determine gdt type from T.
		if(m_band->RasterIO(GF_Read, col, row, cols, rows, block, cols, rows, GDT_Float32, 0, 0) != CE_None)
			throw "Error reading block.";
	}

	/**
	 * Return the resolution. If the x and y resolution are different,
	 * use resolutionX and resolutionY.
	 */
	double resolution() {
		return m_trans[1];
	}

	/**
	 * Get the x resolution.
	 */
	double resolutionX() {
		return m_trans[1];
	}

	/**
	 * Get the y resolution.
	 */
	double resolutionY() {
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

	/**
	 * Returns pixel value at the given coordinate.
	 */
	T get(double x, double y) {
		return get(toCol(x), toRow(y));
	}

	/**
	 * Returns the pixel value at the give row/column.
	 */
	T get(int col, int row) {
		loadBlock(col, row);
		T v = m_block[(row % m_bh) * m_bw + (col % m_bw)];
		return v;
	}

	/**
	 * Sets the pixel value at the given row/column.
	 */
	void set(int col, int row, T v) {
		if(!m_writable) return;
		loadBlock(col, row);
		m_block[(row % m_bh) * m_bw + (col % m_bw)] = v;
		m_dirty = true;
	}

	/**
	 * Sets the pixel value at the given coordinate.
	 */
	void set(double x, double y, T v) {
		set(toCol(x), toRow(y), v);
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
