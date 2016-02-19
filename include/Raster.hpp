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
class Raster {
private:
	int r_cols, r_rows;			// Raster cols/rows
	int r_bcols, r_brows;		// Block cols/rows -- not the number of cols/rows in a block
	int r_curcol, r_currow;		// The current block column
	int r_bandn;				// The band number
	int r_bw, r_bh;				// Block width/height in pixels
	bool r_writable;			// True if the raster is writable
	bool r_dirty;
	T r_nodata;					// Nodata value.
	T *r_block;					// Block storage
	GDALDataset *r_ds;			// GDAL dataset
	GDALRasterBand *r_band;		// GDAL band
	double r_trans[6];			// Raster transform

public:
	Raster() :
		r_cols(-1), r_rows(-1),
		r_bcols(-1), r_brows(-1),
		r_curcol(-1), r_currow(-1),
		r_bandn(1),
		r_bw(-1), r_bh(-1),
		r_writable(false), r_dirty(false) {

		r_ds = nullptr;
		r_band = nullptr;
		r_block = nullptr;

	}
	Raster(std::string &filename, double minx, double miny, double maxx, double maxy,
			double resolution, double nodata, const char *proj = NULL) : Raster() {
		GDALAllRegister();
		int width = (int) ((maxx - minx) / resolution) + 1;
		int height = (int) ((maxy - miny) / resolution) + 1;
		r_ds = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(filename.c_str(),
				width, height, 1, GDT_Float32, NULL);
		if(r_ds == nullptr)
			throw "Failed to create file.";
		// TODO: Proper type.
		r_trans[0] = minx, r_trans[1] = resolution, r_trans[2] = 0.0,
				r_trans[3] = maxy, r_trans[4] = 0.0, r_trans[5] = -resolution;
		r_ds->SetGeoTransform(r_trans);
		if(proj != NULL)
			r_ds->SetProjection(proj);
		r_rows = r_ds->GetRasterYSize();
		r_cols = r_ds->GetRasterXSize();
		r_band = r_ds->GetRasterBand(1);
		if(r_band == NULL)
			throw "Failed to get band.";
		r_band->GetBlockSize(&r_bw, &r_bh);
		r_band->SetNoDataValue(nodata);
		r_nodata = r_band->GetNoDataValue();
		r_bcols = (r_cols + r_bw - 1) / r_bw;
		r_brows = (r_rows + r_bh - 1) / r_bh;
		r_block = (T *) malloc(sizeof(T) * r_bw * r_bh);
		r_writable = true;
	}
	Raster(std::string &filename, int band = 1, bool writable = false) : Raster() {
		GDALAllRegister();
		r_ds = (GDALDataset *) GDALOpen(filename.c_str(), writable ? GA_Update : GA_ReadOnly);
		if(r_ds == NULL)
			throw "Failed to open raster.";
		r_bandn = band;
		r_ds->GetGeoTransform(r_trans);
		r_band = r_ds->GetRasterBand(band);
		if(r_band == nullptr)
			throw "Failed to get band.";
		r_band->GetBlockSize(&r_bw, &r_bh);
		r_rows = r_ds->GetRasterYSize();
		r_cols = r_ds->GetRasterXSize();
		r_nodata = r_band->GetNoDataValue();
		r_bcols = (r_cols + r_bw - 1) / r_bw;
		r_brows = (r_rows + r_bh - 1) / r_bh;
		r_block = (T *) malloc(sizeof(T) * r_bw * r_bh);
		r_writable = writable;
	}
	double resolution() {
		return r_trans[1];
	}
	double resolutionX() {
		return r_trans[1];
	}
	double resolutionY() {
		return r_trans[5];
	}
	void projection(std::string &proj) const {
		proj.assign(r_ds->GetProjectionRef());
	}
	GDALDataType type() const {
		return r_band->GetRasterDataType();
	}
	double minx() const {
		return toX(0);
	}
	double maxx() const {
		return toX(cols());
	}
	double miny() const {
		return toY(rows());
	}
	double maxy() const {
		return toY(0);
	}
	double width() const {
		return maxx() - minx();
	}
	double height() const {
		return maxy() - miny();
	}
	T nodata() const {
		return r_nodata;
	}
	// Returns the row offset in the block for a given y
	int toBlockRow(double y) const {
		return toRow(y) % r_brows;
	}
	// Returns the row offset in the block for a given y
	int toBlockCol(double x) const {
		return toCol(x) % r_bcols;
	}
	// Returns the width of the block
	int blockWidth() const {
		return r_bw;
	}
	// Returns the width of the block
	int blockHeight() const {
		return r_bh;
	}
	// Returns the total number of columns
	int cols() const {
		return r_cols;
	}
	// Returns the total number of rows
	int rows() const {
		return r_rows;
	}
	// Returns the row for a given y-coordinate
	int toRow(double y) const {
		return (int) ((y - r_trans[3]) / r_trans[5]);
	}
	// Returns the column for a given x-coordinate
	int toCol(double x) const {
		return (int) ((x - r_trans[0]) / r_trans[1]);
	}
	// Returns the x-coordinate for a given column
	double toX(int col) const {
		return (col * r_trans[1]) + r_trans[0];
	}
	// Returns the y-coordinate for a given row
	double toY(int row) const {
		return (row * r_trans[5]) + r_trans[3];
	}
	// Returns true if the pixel is nodata
	bool isNoData(int col, int row) {
		return get(col, row) == r_nodata;
	}
	// Returns true if the pixel is nodata
	bool isNoData(double x, double y) {
		return isNoData(toCol(x), toRow(y));
	}
	bool isValid(int c, int r) {
		return getOrNodata(c, r) != r_nodata;
	}
	bool isValid(double x, double y) {
		return getOrNodata(x, y) != r_nodata;
	}
	T getOrNodata(double x, double y) {
		if(!has(x, y)) {
			return r_nodata;
		} else {
			return get(x, y);
		}
	}
	T getOrNodata(int col, int row) {
		if(!has(col, row)) {
			return r_nodata;
		} else {
			return get(col, row);
		}
	}
	// Returns pixel value at the given coordinate
	T get(double x, double y) {
		return get(toCol(x), toRow(y));
	}
	// Returns the pixel value at the give row/column
	T get(int col, int row) {
		loadBlock(col, row);
		T v = r_block[(row % r_bh) * r_bw + (col % r_bw)];
		return v;
	}
	void loadBlock(int col, int row) {
		if(!has(col, row))
			throw "Row or column out of bounds.";
		int bcol = (int) (col / r_bw);
		int brow = (int) (row / r_bh);
		if(bcol >= r_bcols || bcol < 0 || brow >= r_brows || brow < 0)
			throw "Illegal block column or row.";
		if(bcol != r_curcol || brow != r_currow) {
			flush();
			if(r_band->ReadBlock(bcol, brow, r_block) != CE_None)
				throw "Failed to read block.";
			r_currow = brow;
			r_curcol = bcol;
		}
	}
	// Sets the pixel value at the given row/column
	void set(int col, int row, T v) {
		if(!r_writable) return;
		loadBlock(col, row);
		r_block[(row % r_bh) * r_bw + (col % r_bw)] = v;
		r_dirty = true;
	}
	// Sets the pixel value at the given coordinate
	void set(double x, double y, T v) {
		set(toCol(x), toRow(y), v);
	}
	// Flush the current block to the dataset
	void flush() {
		if(r_writable && r_dirty) {
			if(r_band->WriteBlock(r_curcol, r_currow, r_block) != CE_None)
				throw "Flush error.";
			r_dirty = false;
		}
	}
	// Returns true if the col/row are represented in the dataset
	bool has(int col, int row) const {
		return col >= 0 && col < r_cols && row >= 0 && row < r_rows;
	}
	bool has(double x, double y) const {
		return has(toCol(x), toRow(y));
	}
	~Raster() {
		flush();
		r_band = NULL;
		r_ds = NULL;
		std::cerr << "Flushed" << std::endl;
	}

};


#endif /* INCLUDE_RASTER_HPP_ */
