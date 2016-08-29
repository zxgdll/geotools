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
#include <string>

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
			Cell(int col, int row);
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
			TargetOperator(T match);

			bool fill(T value) const;
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
			void gaussianWeights(double *weights, int size, double sigma);

		public:
			/**
			 * Return the number of rows in the dataset.
			*/
			virtual int rows() const =0;

			/**
			 * Return the number of columns in the dataset.
			 */
			virtual int cols() const =0;
			
			/**
			 * Return the number of cells in the dataset.
			 */
			virtual size_t size() const =0;

			/**
			 * Fill the entire dataset with the given value.
			 */
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
			 * Return a reference to the value held at the given index in the grid.
			 * Not const because the get operation might imply (e.g.) a buffering operation in the subclass.
			 */
			virtual T &get(size_t idx) =0;

			/**
			 * Return a reference to the value held at the given column and row.
			 * Not const because the get operation might imply (e.g.) a buffering operation in the subclass.
			 */
			virtual T &get(int col, int row) =0;

			/**
			 * Set the value held at  the given index in the grid.
			 */	
			virtual void set(size_t idx, const T value) =0;

			/**
			 * Set the value held at  the given column and row.
			 */	
			virtual void set(int col, int row, const T value) =0;

			/**
			 * Get or set the value held at the given index in the grid.
			 * Not const because the get operation might imply (e.g.) a buffering operation in the subclass.
			 */
			virtual T &operator[](size_t idx) =0;
			
			/* 
			 * Return true if the dataset contains the given element.
			 */
			virtual bool has(int col, int row) const =0;
			
			/* 
			 * Return true if the dataset contains the given element.
			 */
			virtual bool has(size_t idx) const =0;

			/**
			 * Returns trueif the dataset contains the given element and it is nodata.
			 */
			virtual bool isNoData(int col, int row) =0;
			
			/**
			 * Returns trueif the dataset contains the given element and it is nodata.
			 */
			virtual bool isNoData(size_t idx) =0;

			/**
			 * Returns true if the grid is a perfect square.
			 */
			virtual bool isSquare() const =0;

			/**
			 * Returns the nodata value.
			 */
			virtual T nodata() const =0;

			/**
			 * Sets the nodata value.
			 */
			virtual void nodata(T nodata) =0;

			/**
			 * Read data into Grid instance. Will attempt to read a region of the same size
			 * as the given block.
			 * col    - The column in the data source to read from.
			 * row    - The row in the data source to read from.
			 * block  - The block to read into.
			 * dstCol - The column in the destination block to write to.
			 * dstRow - The row in the destination block to write to. 
			 */
			virtual void readBlock(int col, int row, Grid<T> &block, int dstCol = 0, int dstRow = 0) =0;

			/**
			 * Read data into Grid instance. Will attempt to read a region of the same size
			 * as the given block.
			 */
			virtual void readBlock(Grid<T> &block) =0;
			
			/**
			 * Write data from Grid instance. Will attempt to write a region of the same size
			 * as the given block.
			 * col    - The column in the data source to write to..
			 * row    - The row in the data source to write to.
			 * block  - The block to write from..
			 * srcCol - The column in the source block to read from.
			 * srcRow - The row in the source block to read from. 
			 */
			virtual void writeBlock(int col, int row, Grid<T> &block, int dstCol = 0, int dstRow = 0) =0;

			/**
			 * Write data from Grid instance. Will attempt to write a region of the same size
			 * as the given block.
			 */
			virtual void writeBlock(Grid<T> &block) =0;


			/**
			 * Computes descriptive statistics for the values in the grid.
			 * A nodata value must be set.
			 */
			void computeStats();

			/**
			 * Return the maximum value in the raster.
			 */
			T max();

			/**
			 * Return the minimum value in the raster.
			 */
			T min();

			/**
			 * Return the mean value in the raster.
			 */
			T mean();

			/**
			 * Return the standard deviation of  values in the raster.
			 */	
			T stddev();

			/**
			 * Return the variance of values in the raster.
			 */
			T variance();

			/**
			 * Fill the grid, beginning with the target cell, where any contiguous cell
			 * satisfies the given FillOperator. The other grid is actually filled,
			 * and the present grid is unchanged *unless* the present grid is passed
			 * as other.
			 */
			template <class U>
			std::vector<int> floodFill(int col, int row, FillOperator<T> &op, Grid<U> &other, U fill);

			/**
			 * Begin flood fill at the given cell; fill cells equal to the target value.
			 */
			std::vector<int> floodFill(int col, int row, T target, T fill);

			/**
			 * Begin flood fill at the given cell; fill cells that satisfy the operator.
			 */
			std::vector<int> floodFill(int col, int row, FillOperator<T> &op, T fill);

			/**
			 * Smooth the raster and write the smoothed version to the output raster.
			 */
			void smooth(Grid<T> &smoothed, double sigma, int size);

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
			void checkInit() const;

		public:
			MemRaster();

			MemRaster(int cols, int rows);

			template <class U>
			MemRaster(Grid<U> &tpl);

			~MemRaster();

			/**
			 * A pointer to a function that can deallocate the grid
			 * items. Not necessary if a primitive type is used (the usual case.)
			 */
			void setDeallocator(void (*item_dealloc)(T));

			/**
			 * Return a pointer to the allocated memory.
			 */
			T *grid();

			bool hasGrid() const;

			/**
			 * Cas a MemRaster to some other type.
			 */
			template <class U>
			void convert(MemRaster<U> &g);

			int rows() const;

			int cols() const;

			size_t size() const;

			template <class U>
			void init(Grid<U> &tpl);

			/**
			 * Initialize with the given number of cols and rows.
			 * (Re)allocates memory for the internal grid.
			 */
			void init(int cols, int rows);

			/**
			 * Fill the grid with the given value.
			 */
			void fill(const T value);

			/**
			 * Return a reference to the value held at
			 * the given index in the grid.
			 */
			T &get(size_t idx);

			T &get(int col, int row);

			bool isNoData(int col, int row);

			bool isNoData(size_t idx);

			void set(int col, int row, const T value);

			void set(size_t idx, const T value);

			bool has(int col, int row) const;

			bool has(size_t idx) const;

			/**
			 * Return the element at the given index.
			 */
			T &operator[](size_t idx);

			bool isSquare() const;

			/**
			 * Convert the grid to matrix.
			 */
			void toMatrix(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx);

			/**
			 * Initialize the grid from a matrix.
			 */
			void fromMatrix(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx);

			T nodata() const;

			void nodata(T nodata);

			void readBlock(int col, int row, Grid<T> &block, int dstCol = 0, int dstRow = 0);

			void writeBlock(int col, int row, Grid<T> &block, int srcCol = 0, int srcRow = 0);

			void writeBlock(Grid<T> &block);

			void readBlock(Grid<T> &block);

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

			/**
			 * Flush the block at the given index to disk, if it is dirty.
			 */
			void flushBlock(size_t idx);

			/**
			 * Convert the row and column indices to an index.
			 */
			size_t toIdx(int col, int row);

			/**
			 * Convert the index to a column index.
			 */
			int toCol(size_t idx);

			/**
			 * Convert the index to a row index.
			 */
			int toRow(size_t idx);

			/**
			 * Remove the oldest block from the cache and return a reference
			 * to its memory. Flush first if dirty.
			 */
			T* freeOldest();

			/**
			 * Remove a block from the cache and return a reference to its
			 * memory for reuse. If there are more than the limit of blocks
			 * keep removing until the cache is full-1.
			 * Will flush dirty blocks.
			 */
			T* freeOne();

		public:
			BlockCache();

			/**
			 * The number of columns in a single block.
			 */
			int blockWidth();

			/**
			 * The number of rows in a single block.
			 */
			int blockHeight();

			/**
			 * Convert the raster row and column into a row and
			 * column index for a single block.
			 */
			size_t toBlockIdx(int col, int row);
			
			/**
			 * Set the raster band for this cache to manage.
			 */
			void setRasterBand(GDALRasterBand *band);

			/**
			 * Return true if the cache is managing a block that contains the
			 * given column and row.
			 */
			bool hasBlock(int col, int row);

			/**
			 * Return true if the given index is valid for this cache.
			 */
			bool hasBlock(size_t idx);

			/**
			 * Set the number of blocks managed by the cache.
			 */
			void setSize(size_t size);

			/**
			 * Return the number of blocks managed by the cache.
			 */
			size_t getSize();

			/**
			 * Return a pointer to the block containing the given column and row.
			 * If the block is to be written, the forWrite argument flags it as 
			 * dirty.
			 */
			T* getBlock(int col, int row, bool forWrite);

			/**
			 * Flush all blocks to disk.
			 */
			void flush();

			/**
			 * Free all blocks.
			 */
			void close();

			~BlockCache();

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
			void loadBlock(int col, int row, bool forWrite);

			/**
			 * Get the GDAL type for the given c++ type.
			 */
			GDALDataType getType(double v);

			GDALDataType getType(float v);

			GDALDataType getType(unsigned int v);

			GDALDataType getType(unsigned long v);
			
			GDALDataType getType(int v);

			GDALDataType getType(unsigned short v);

			GDALDataType getType(short v);

			GDALDataType getType(char v);

			/**
			 * Get the GDAL type for the current raster.
			 */
			GDALDataType getType();

			/**
			 * Return the default nodata for the raster's type.
			 */
			T getDefaultNodata();

		public:

			/**
			 * Basic constructor.
			 */
			Raster();

			/**
			 * Create a new raster for writing with a template of a different type.
			 */
			template <class U>
			Raster(const std::string &filename, int band, const Raster<U> &tpl);

			/**
			 * Create a new raster for writing with a template of a different type.
			 */
			//Raster(const std::string &filename, int band, const Raster<T> &tpl);

			/**
			 * Build a new raster with the given filename, bounds, resolution, nodata and projection.
			 */
			Raster(const std::string &filename, int band, double minx, double miny, double maxx, double maxy,
					double resolutionX, double resolutionY, double nodata, const std::string &proj);

			/**
			 * Build a new raster with the given filename, bounds, resolution, nodata and projection.
			 */
			Raster(const std::string &filename, int band, Bounds &bounds, double resolutionX, double resolutionY, double nodata, int crs);
			
			/**
			 * Build a new raster with the given filename, bounds, resolution, nodata and SRID.
			 */
			Raster(const std::string &filename, int band, double minx, double miny, double maxx, double maxy,
					double resolutionX, double resolutionY, double nodata, int crs);

			/**
			 * Open the given raster and load the given band. Set the writable argument to true
			 * to enable writing.
			 */
			Raster(const std::string &filename, int band = 1, bool writable = false);

			/**
			 * Initialize the raster using the given file (which may not exist) using
			 * another raster as a template. The raster pixel types need not be the same.
			 */
			template <class U>
			void init(const std::string &filename, int band, const Raster<U> &tpl);

			void init(const std::string &filename, int band, const Bounds &bounds, double resolutionX, double resolutionY,
				double nodata, const std::string &proj);

			/**
			 * Initializes the raster with the given filename, bounds, resolution, nodata and projection.
			 */
			void init(const std::string &filename, int band, double minx, double miny, double maxx, double maxy,
					double resolutionX, double resolutionY, double nodata, const std::string &proj);

			/**
			 * Initializes a Raster from the existing file.
			 */
			void init(const std::string &filename, int band = 1, bool writable = false);

			/**
			 * Set the number of blocks in the cache.
			 */
			void setCacheSize(size_t size);

			/**
			 * Return the filename for this raster.
			 */
			std::string filename() const;

			/**
			 * Return the number of bands in the raster.
			 */
			int bandCount() const;

			/**
			 * Converts a numerical (EPSG) crs code to a projection string.
			 */
			std::string epsg2ProjText(int crs) const;

			/**
			 * Returns true if the raster is initialized.
			 */
			bool inited() const;

			void fill(T value);

			void readBlock(int col, int row, Grid<T> &grd, int dstCol = 0, int dstRow = 0);

			void writeBlock(int col, int row, Grid<T> &grd, int srcCol = 0, int srcRow = 0);

			void writeBlock(Grid<T> &block);
			
			void readBlock(Grid<T> &block);
			
			/**
			 * Get the x resolution.
			 */
			double resolutionX() const;

			/**
			 * Get the y resolution.
			 */
			double resolutionY() const;

			/**
			 * Write the projection data to the given string object.
			 */
			void projection(std::string &proj) const;

			/**
			 * Return the GDAL datatype of the raster.
			 */
			GDALDataType type() const;

			/**
			 * The minimum bounding x.
			 */
			double minx() const;

			/**
			 * The maximum bounding x.
			 */
			double maxx() const;

			/**
			 * The minimum bounding y.
			 */
			double miny() const;

			/**
			 * The maximum bounding y.
			 */
			double maxy() const;

			/**
			 * The width of the raster in map units.
			 */
			double width() const;

			/**
			 * The height of the raster in map units.
			 */
			double height() const;

			T nodata() const;

			void nodata(T nodata);

			int cols() const;

			int rows() const;

			/**
			 * Returns the row for a given y-coordinate.
			 */
			int toRow(double y) const;

			/**
			 * Returns the column for a given x-coordinate.
			 */
			int toCol(double x) const;

			/**
			 * Returns the x-coordinate for a given column.
			 */
			double toX(int col) const;

			/**
			 * Returns the y-coordinate for a given row.
			 */
			double toY(int row) const;

			size_t size() const;

			/**
			 * Returns true if the pixel is nodata.
			 */
			bool isNoData(int col, int row);

			/**
			 * Returns true if the pixel is nodata.
			 */
			bool isNoData(size_t idx);

			/**
			 * Returns true if the pixel is nodata.
			 */
			bool isNoData(double x, double y);

			/**
			 * Returns true if the pixel exists and is not nodata.
			 */
			bool isValid(int c, int r);

			/**
			 * Returns true if the pixel exists and is not nodata.
			 */
			bool isValid(double x, double y);

			/**
			 * Gets the pixel value or nodata if the pixel doesn't exist.
			 */
			T getOrNodata(double x, double y);

			/**
			 * Gets the pixel value or nodata if the pixel doesn't exist.
			 */
			T getOrNodata(int col, int row);

			T *grid();
			
			bool hasGrid() const;

			T &get(double x, double y);

			T &get(int col, int row);

			T &get(size_t idx);

			T &operator[](size_t idx);

			void set(int col, int row, T v);

			void set(size_t idx, T v);

			void set(double x, double y, T v);

			bool isSquare() const;

			bool has(int col, int row) const;

			bool has(double x, double y) const;

			bool has(size_t idx) const;

			/**
			 * Flush the current block to the dataset.
			 */
			void flush();

			~Raster();

		};

	} // raster

} // geotools


#endif
