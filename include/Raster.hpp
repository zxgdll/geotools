/*
 * Raster.hpp
 *
 *  Created on: Jan 21, 2016
 *      Author: rob
 */

#ifndef INCLUDE_RASTER_HPP_
#define INCLUDE_RASTER_HPP_

template <class T>
class Raster {
private:
	int r_cols, r_rows;			// Raster cols/rows
	int r_bcols, r_brows;		// Block cols/rows -- not the number of cols/rows in a block
	int r_curcol, r_currow;		// The current block column
	int r_bandn;				// The band number
	int r_bw, r_bh;				// Block width/height in pixels
	bool r_writable;			// True if the raster is writable
	T r_nodata;					// Nodata value.
	GDALDataset *r_ds;			// GDAL dataset
	GDALRasterBand *r_band;		// GDAL band
	double r_trans[6];			// Raster transform
	T *r_block;					// Block storage
	bool r_close;					// If true, close handle on destruct.

public:
	Raster(std::string &filename, int band = 1, bool writable = false) {
		GDALAllRegister();
		r_ds = (GDALDataset *) GDALOpen(filename.c_str(), writable ? GA_Update : GA_ReadOnly);
		if(r_ds == NULL)
			throw "Failed to open raster.";
		r_bandn = band;
		r_curcol = -1;
		r_currow = -1;
		r_ds->GetGeoTransform(r_trans);
		r_band = r_ds->GetRasterBand(band);
		r_band->GetBlockSize(&r_bw, &r_bh);
		r_rows = r_ds->GetRasterYSize();
		r_cols = r_ds->GetRasterXSize();
		r_nodata = r_band->GetNoDataValue();
		r_bcols = r_cols / r_bw;
		r_brows = r_rows / r_bh;
		r_block = (T *) malloc(sizeof(T) * r_bw * r_bh);
		r_writable = writable;
		r_close = true;
	}
	Raster(GDALDataset *ds, int band = 1, bool writable = false) {
		r_ds = ds;
		r_bandn = band;
		r_curcol = -1;
		r_currow = -1;
		ds->GetGeoTransform(r_trans);
		r_band = ds->GetRasterBand(band);
		r_band->GetBlockSize(&r_bw, &r_bh);
		r_rows = ds->GetRasterYSize();
		r_cols = ds->GetRasterXSize();
		r_nodata = r_band->GetNoDataValue();
		r_bcols = r_cols / r_bw;
		r_brows = r_rows / r_bh;
		r_block = (T *) malloc(sizeof(T) * r_bw * r_bh);
		r_writable = writable;
		r_close = false;
	}
	Raster<T> copy(std::string filename, int band = 1, bool writable = false) {
		GDALDataset *outDS = r_ds->GetDriver()->CreateCopy(filename.c_str(), r_ds, band, NULL, NULL, NULL);
		if(outDS == NULL)
			throw "Failed to copy file.";
		GDALClose(outDS);
		Raster<T> cp(filename, band, writable);
		T nodata = cp.nodata();
		for(int r = 0; r < cp.rows(); ++r) {
			for(int c = 0; c < cp.cols(); ++c) {
				cp.set(c, r, nodata);
			}
		}
		cp.flush();
		return cp;
	}
	T nodata() {
		return r_nodata;
	}
	// Returns the row offset in the block for a given y
	int toBlockRow(float y) {
		return toRow(y) % r_brows;
	}
	// Returns the row offset in the block for a given y
	int toBlockCol(float x) {
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
	int toRow(float y) const {
		return (int) ((y - r_trans[3]) / r_trans[5]);
	}
	// Returns the column for a given x-coordinate
	int toCol(float x) const {
		return (int) ((x - r_trans[0]) / r_trans[1]);
	}
	// Returns the x-coordinate for a given column
	float toX(int col) const {
		return (col * r_trans[1]) + r_trans[0];
	}
	// Returns the y-coordinate for a given row
	float toY(int row) const {
		return (row * r_trans[5]) + r_trans[3];
	}
	// Returns true if the pixel is nodata
	bool isNoData(int col, int row) {
		return get(col, row) == r_nodata;
	}
	// Returns true if the pixel is nodata
	bool isNoData(float x, float y) {
		return isNoData(toCol(x), toRow(y));
	}
	// Returns pixel value at the given coordinate
	T get(float x, float y) {
		return get(toCol(x), toRow(y));
	}
	// Returns the pixel value at the give row/column
	T get(int col, int row) {
		if(!has(col, row))
			throw "Row or column out of bounds.";
		int bcol = (int) (col / r_bw);
		int brow = (int) (row / r_bh);
		if(bcol != r_curcol || brow != r_currow) {
			//std::cerr << "New block " << bcol << ", " << brow << std::endl;
			flush();
			CPLErr ret = r_band->ReadBlock(bcol, brow, r_block);
			if(ret != 0)
				throw "Failed to read block.";
			r_currow = brow;
			r_curcol = bcol;
		}
		return r_block[(row % r_bh) * r_bw + (col % r_bw)];
	}
	// Sets the pixel value at the given row/column
	void set(int col, int row, T v) {
		if(!r_writable) return;
		get(col, row);
		r_block[(row % r_bh) * r_bw + (col % r_bw)] = v;
	}
	// Sets the pixel value at the given coordinate
	void set(float x, float y, T v) {
		set(toCol(x), toRow(y), v);
	}
	// Flush the current block to the dataset
	void flush() {
		if(r_writable && r_curcol > -1 && r_currow > -1) {
			GDALRasterBlock *b = r_band->GetLockedBlockRef(r_curcol, r_currow);
			CPLErr ret = r_band->WriteBlock(r_curcol, r_currow, r_block);
			if(ret != 0)
				throw "Flush error.";
			b->DropLock();
		}
	}
	// Returns true if the col/row are represented in the dataset
	bool has(int col, int row) const {
		return col >= 0 || col < r_cols || row >= 0 || row < r_rows;
	}
	void close() {
		flush();
		free(r_block);
		r_ds = 0;
		r_band = 0;
		if(r_close && r_ds)
			GDALClose(r_ds);
	}
	~Raster() {
		close();
	}

};


#endif /* INCLUDE_RASTER_HPP_ */
