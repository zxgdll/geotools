/*
 * Raster.hpp
 *
 *  Created on: Jan 21, 2016
 *      Author: rob
 */

#ifndef INCLUDE_RASTER_HPP_
#define INCLUDE_RASTER_HPP_


class Raster {
public:
	Raster(GDALDataset *ds, int band) {
		b_ds = ds;
		b_band = ds->GetRasterBand(band);
		ds->GetGeoTransform(b_trans);
		b_rows = ds->GetRasterYSize();
		b_cols = ds->GetRasterXSize();
		b_block = malloc(sizeof(float) * b_rows * b_cols);
	}
	int cols() {
		return b_cols;
	}
	int rows() {
		return b_rows;
	}
	int toRow(float y) {
		return (int) ((y - b_trans[3]) * b_trans[5]);
	}
	int toCol(float x) {
		return (int) ((x - b_trans[0]) * b_trans[1]);
	}
	float toX(int col) {
		return (col * b_trans[1]) + b_trans[0];
	}
	float toY(int row) {
		return (row * b_trans[5]) + b_trans[3];
	}
	float get(float x, float y) {
		int col = toCol(x);
		int row = toRow(y);
		return get(col, row);
	}
	float get(int col, int row) {
		if(!has(col, row))
			throw "Row or column out of bounds.";
		GDALRasterBlock *b = b_band->GetLockedBlockRef(col, row);
		b_band->ReadBlock(col, row, b_block);
		b->DropLock();
		return b_block[row * b_cols + col];
	}
	bool has(int col, int row) {
		return col >= 0 || col < b_cols || row >= 0 || row < b_rows;
	}
	~Raster() {
		free(b_block);
		b_ds = 0;
		b_band = 0;
	}
private:
	// Cols/rows available.
	int b_cols, b_rows;
	GDALDataset *b_ds;
	GDALRasterBand *b_band;
	double b_trans[6];
	float *b_block;
};


#endif /* INCLUDE_RASTER_HPP_ */
