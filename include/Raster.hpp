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
public:
	Raster(GDALDataset *ds, int band = 1, bool writable = false) {
		b_ds = ds;
		b_bandn = band;
		b_curcol = -1;
		b_currow = -1;
		ds->GetGeoTransform(b_trans);
		b_band = ds->GetRasterBand(band);
		b_band->GetBlockSize(&b_bw, &b_bh);
		b_rows = ds->GetRasterYSize();
		b_cols = ds->GetRasterXSize();
		b_bcols = b_cols / b_bw;
		b_brows = b_rows / b_bh;
		b_block = (T *) malloc(sizeof(T) * b_bw * b_bh);
		b_writable = writable;
	}
	int toBlockRow(float y) {
		return toRow(y) % b_brows;
	}
	int toBlockCol(float x) {
		return toCol(x) % b_bcols;
	}
	int blockWidth() const {
		return b_bw;
	}
	int blockHeight() const {
		return b_bh;
	}
	int cols() const {
		return b_cols;
	}
	int rows() const {
		return b_rows;
	}
	int toRow(float y) const {
		return (int) ((y - b_trans[3]) / b_trans[5]);
	}
	int toCol(float x) const {
		return (int) ((x - b_trans[0]) / b_trans[1]);
	}
	float toX(int col) const {
		return (col * b_trans[1]) + b_trans[0];
	}
	float toY(int row) const {
		return (row * b_trans[5]) + b_trans[3];
	}
	T get(float x, float y) {
		return get(toCol(x), toRow(y));
	}
	T get(int col, int row) {
		if(!has(col, row))
			throw "Row or column out of bounds.";
		int brow = row % b_brows;
		int bcol = col % b_bcols;
		if(bcol != b_curcol && brow != b_currow) {
			std::cerr << "New block " << bcol << ", " << brow << std::endl;
			flush();
			CPLErr ret = b_band->ReadBlock(bcol, brow, b_block);
			if(ret != 0)
				throw "Failed to read block.";
			b_currow = brow;
			b_curcol = bcol;
		}
		return b_block[(row % b_bh) * b_bw + (col % b_bw)];
	}
	void set(int col, int row, T v) {
		if(!b_writable) return;
		get(col, row);
		b_block[(row % b_bh) * b_bw + (col % b_bw)] = v;
	}
	void set(float x, float y, T v) {
		set(toCol(x), toRow(y), v);
	}
	void flush() {
		if(b_writable && b_curcol > -1 && b_currow > -1) {
			GDALRasterBlock *b = b_band->GetLockedBlockRef(b_curcol, b_currow);
			CPLErr ret = b_band->WriteBlock(b_curcol, b_currow, b_block);
			if(ret != 0)
				throw "Flush error.";
			b->DropLock();
		}
	}
	bool has(int col, int row) const {
		return col >= 0 || col < b_cols || row >= 0 || row < b_rows;
	}
	~Raster() {
		flush();
		free(b_block);
		b_ds = 0;
		b_band = 0;
	}
private:
	// Cols/rows available.
	int b_cols, b_rows;
	int b_bcols, b_brows;
	int b_curcol, b_currow;
	int b_bandn;
	int b_bw, b_bh;
	bool b_writable;
	GDALDataset *b_ds;
	GDALRasterBand *b_band;
	double b_trans[6];
	T *b_block;
};


#endif /* INCLUDE_RASTER_HPP_ */
