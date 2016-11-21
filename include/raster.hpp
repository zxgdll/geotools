/**
 * Raster.hpp
 *
 *  Created on: Jan 21, 2016
 *  Author: rob
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
#include <memory>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <eigen3/Eigen/Core>

#include "geotools.h"
#include "util.hpp"

using namespace geotools::util;

namespace geotools {

    namespace raster {

        // Simple class to represent a single grid cell.

        class DLL_EXPORT Cell {
        public:
            int32_t col;
            int32_t row;
            Cell(int32_t col, int32_t row);
        };

        // Used by Grid::floodFill to determine whether
        // a pixel should be filled.

        template <class T>
        class DLL_EXPORT FillOperator {
        public:
            virtual bool fill(T value) const = 0;
        };

        template <class T>
        class DLL_EXPORT TargetOperator : public FillOperator<T> {
        private:
            T m_match;
        public:
            TargetOperator(T match);

            bool fill(T value) const;
        };

        // Abstract class for grids (rasters).

        template <class T>
        class DLL_EXPORT Grid {
        protected:
            T m_min;
            T m_max;
            T m_mean;
            T m_stddev;
            T m_variance;
            T m_sum;
            int32_t m_count;
            bool m_stats;

            // Compute the table of Gaussian weights given the size of the table
            // and the std. deviation.
            void gaussianWeights(double *weights, int32_t size, double sigma);

        public:
            Grid();

            // Return the number of rows in the dataset.
            virtual int32_t rows() const = 0;

            // Return the number of columns in the dataset.
            virtual int32_t cols() const = 0;

            // Return the number of cells in the dataset.
            virtual size_t size() const = 0;

            // Fill the entire dataset with the given value.
            virtual void fill(const T value) = 0;

            // Return a pointer to an in-memory grid of the raster data.
            // Throw an exception if this is not possible.
            virtual T *grid() = 0;

            // Returns true if this class has a complete, in-memory
            // grid that can be manipulated.
            virtual bool hasGrid() const = 0;

            // Return a reference to the value held at the given index in the grid.
            // Not const because the get operation might imply (e.g.) a buffering 
            // operation in the subclass.
            virtual T get(size_t idx) = 0;

            // Return a reference to the value held at the given column and row.
            // Not const because the get operation might imply (e.g.) a buffering 
            // operation in the subclass.
            virtual T get(int32_t col, int32_t row) = 0;

            // Set the value held at  the given index in the grid.
            virtual void set(size_t idx, const T value) = 0;

            // Set the value held at  the given column and row.
            virtual void set(int32_t col, int32_t row, const T value) = 0;

            // Get or set the value held at the given index in the grid.
            // Not const because the get operation might imply (e.g.) a buffering 
            // operation in the subclass.
            virtual T operator[](size_t idx) = 0;

            // Return true if the dataset contains the given element.
            virtual bool has(int32_t col, int32_t row) const = 0;

            // Return true if the dataset contains the given element.
            virtual bool has(size_t idx) const = 0;

            // Returns trueif the dataset contains the given element and it is nodata.
            virtual bool isNoData(int32_t col, int32_t row) = 0;

            // Returns trueif the dataset contains the given element and it is nodata.
            virtual bool isNoData(size_t idx) = 0;

            // Returns true if the grid is a perfect square.
            virtual bool isSquare() const = 0;

            // Returns the nodata value.
            virtual T nodata() const = 0;

            // Sets the nodata value.
            virtual void nodata(const T nodata) = 0;

            // Read data into Grid instance. Will attempt to read a region of 
            // the same size as the given block.
            // col    - The column in the data source to read from.
            // row    - The row in the data source to read from.
            // block  - The block to read into.
            // dstCol - The column in the destination block to write to.
            // dstRow - The row in the destination block to write to. 
            // xcols  - The max number of cols to write.
            // xrows  - The max number of rows to write.
            virtual void readBlock(int32_t col, int32_t row, Grid<T> &block, 
                int32_t dstCol = 0, int32_t dstRow = 0, int32_t xcols = 0, 
                int32_t xrows = 0) = 0;

            // Read data into Grid instance. Will attempt to read a region of the same size
            // as the given block.
            virtual void readBlock(Grid<T> &block) = 0;

            // Write data from Grid instance. Will attempt to write a region of the same size
            // as the given block.
            // col    - The column in the data source to write to..
            // row    - The row in the data source to write to.
            // block  - The block to write from..
            // srcCol - The column in the source block to read from.
            // srcRow - The row in the source block to read from. 
            // xcols  - The max number of cols to write.
            // xrows  - The max number of rows to write.
            virtual void writeBlock(int32_t col, int32_t row, Grid<T> &block, 
                int32_t dstCol = 0, int32_t dstRow = 0, int32_t xcols = 0, 
                int32_t xrows = 0) = 0;

            // Write data from Grid instance. Will attempt to write a region of 
            // the same size as the given block.
            virtual void writeBlock(Grid<T> &block) = 0;


            // Computes descriptive statistics for the values in the grid.
            // A nodata value must be set.
            void computeStats();

            // Normalize the grid so that one standard deviation is +-1.
            void normalize();

            // Return the maximum value in the raster.
            T max();

            // Return the minimum value in the raster.
            T min();

            // Return the mean value in the raster.
            T mean();

            // Return the standard deviation of  values in the raster.
            T stddev();

            // Return the variance of values in the raster.
            T variance();

            // Convert a Grid to some other type.
            template <class U>
            void convert(Grid<U> &g) {
                for (size_t i = 0; i < size(); ++i)
                    g.set(i, (U) get(i));
            }

            // Fill the grid, beginning with the target cell, where any contiguous cell
            // satisfies the given FillOperator. The other grid is actually filled,
            // and the present grid is unchanged *unless* the present grid is passed
            // as other.
            // TODO: Moved here to allow compilation of different type combinations.
            template <class U>
            std::vector<int> floodFill(int32_t col, int32_t row, 
                FillOperator<T> &op, Grid<U> &other, U fill) {

                int minc = cols() + 1;
                int minr = rows() + 1;
                int maxc = -1;
                int maxr = -1;
                int area = 0;

                std::queue<std::unique_ptr<Cell> > q;
                q.push(std::unique_ptr<Cell>(new Cell(col, row)));

                std::vector<bool> visited(size(), false); // Tracks visited pixels.

                while (q.size()) {

                    std::unique_ptr<Cell> cel = std::move(q.front());
                    row = cel->row;
                    col = cel->col;
                    q.pop();

                    size_t idx = (size_t) row * cols() + col;

                    if (!visited[idx] && op.fill(get(col, row))) {

                        minc = g_min(col, minc);
                        maxc = g_max(col, maxc);
                        minr = g_min(row, minr);
                        maxr = g_max(row, maxr);
                        ++area;
                        other.set(col, row, fill);
                        visited[idx] = true;

                        if (row > 0)
                            q.push(std::unique_ptr<Cell>(new Cell(col, row - 1)));
                        if (row < rows() - 1)
                            q.push(std::unique_ptr<Cell>(new Cell(col, row + 1)));

                        for (int32_t c = col - 1; c >= 0; --c) {
                            idx = (size_t) row * cols() + c;
                            if (!visited[idx] && op.fill(get(c, row))) {
                                minc = g_min(c, minc);
                                ++area;
                                other.set(c, row, fill);
                                visited[idx] = true;
                                if (row > 0)
                                    q.push(std::unique_ptr<Cell>(new Cell(c, row - 1)));
                                if (row < rows() - 1)
                                    q.push(std::unique_ptr<Cell>(new Cell(c, row + 1)));
                            } else {
                                break;
                            }
                        }
                        for (int32_t c = col + 1; c < cols(); ++c) {
                            idx = (size_t) row * cols() + c;
                            if (!visited[idx] && op.fill(get(c, row))) {
                                maxc = g_max(c, maxc);
                                ++area;
                                other.set(c, row, fill);
                                visited[idx] = true;
                                if (row > 0)
                                    q.push(std::unique_ptr<Cell>(new Cell(c, row - 1)));
                                if (row < rows() - 1)
                                    q.push(std::unique_ptr<Cell>(new Cell(c, row + 1)));
                            } else {
                                break;
                            }
                        }
                    }
                }
                std::vector<int> ret(5);
                ret.push_back(minc);
                ret.push_back(minr);
                ret.push_back(maxc);
                ret.push_back(maxr);
                ret.push_back(area);
                return ret;
            }

            // Begin flood fill at the given cell; fill cells equal to the target value.
            std::vector<int> floodFill(int32_t col, int32_t row, T target, T fill);

            // Begin flood fill at the given cell; fill cells that satisfy the operator.
            std::vector<int> floodFill(int32_t col, int32_t row, FillOperator<T> &op, T fill);

            // Smooth the raster and write the smoothed version to the output raster.
            // Callback is an optional function reference with a single float
            // between 0 and 1, for status tracking.
            void smooth(Grid<T> &smoothed, double sigma, int32_t size, geotools::util::Callbacks *status = nullptr);

            // The radius is given with cells as the unit, but
            // can be rational. When determining which cells to
            // include in the calculation, any cell which partially
            // falls in the radius will be included.
            void voidFillIDW(double radius, int32_t count = 4, double exp = 2.0);

        };

        // A convenience class for managing a grid of values.
        // Handles allocation and deallocation of memory.

        template <class T>
        class DLL_EXPORT MemRaster : public Grid<T> {
        private:
            T *m_grid;
            int32_t m_cols;
            int32_t m_rows;
            void (*m_item_dealloc)(T);
            T m_nodata;
            bool m_mmapped;
            size_t m_size;
            std::unique_ptr<geotools::util::MappedFile> m_mappedFile;
            std::unique_ptr<boost::interprocess::mapped_region> m_region;
            std::unique_ptr<boost::interprocess::file_mapping> m_mapping;

            // Checks if the grid has been initialized. Throws exception otherwise.
            void checkInit() const;

            void freeMem();

        public:
            MemRaster();

            MemRaster(int32_t cols, int32_t rows, bool mapped = false);

            template <class U>
            MemRaster(Grid<U> &tpl, bool mapped = false) : MemRaster() {
                init(tpl.cols(), tpl.rows(), mapped);
            }

            ~MemRaster();

            // A pointer to a function that can deallocate the grid
            // items. Not necessary if a primitive type is used (the usual case.)
            void setDeallocator(void (*item_dealloc)(T));

            // Return a pointer to the allocated memory.
            T *grid();

            bool hasGrid() const;

            int32_t rows() const;

            int32_t cols() const;

            size_t size() const;

            template <class U>
            void init(Grid<U> &tpl, bool mapped = false) {
                init(tpl.cols(), tpl.rows(), mapped);
            }

            // Initialize with the given number of cols and rows.
            // (Re)allocates memory for the internal grid.
            void init(int32_t cols, int32_t rows, bool mapped = false);

            // Fill the grid with the given value.
            void fill(const T value);

            // Return a reference to the value held at
            // the given index in the grid.
            T get(size_t idx);

            T get(int32_t col, int32_t row);

            bool isNoData(int32_t col, int32_t row);

            bool isNoData(size_t idx);

            void set(int32_t col, int32_t row, const T value);

            void set(size_t idx, const T value);

            bool has(int32_t col, int32_t row) const;

            bool has(size_t idx) const;

            // Return the element at the given index.
            T operator[](size_t idx);

            bool isSquare() const;

            // Convert the grid to matrix.
            void toMatrix(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx);

            // Initialize the grid from a matrix.
            void fromMatrix(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx);

            T nodata() const;

            void nodata(const T nodata);

            void readBlock(int32_t col, int32_t row, Grid<T> &block, int32_t dstCol = 0, int32_t dstRow = 0, int32_t xcols = 0, int32_t xrows = 0);

            void writeBlock(int32_t col, int32_t row, Grid<T> &block, int32_t srcCol = 0, int32_t srcRow = 0, int32_t xcols = 0, int32_t xrows = 0);

            void writeBlock(Grid<T> &block);

            void readBlock(Grid<T> &block);

        };

        template <class T>
        class DLL_EXPORT BlockCache {
        private:
            GDALRasterBand *m_band;
            size_t m_size;
            size_t m_time;
            int32_t m_bw;
            int32_t m_bh;
            std::unordered_map<size_t, T*> m_blocks; // idx, block
            std::unordered_map<size_t, size_t> m_idx_time; // idx, time
            std::unordered_map<size_t, bool> m_dirty; // idx, bool
            std::map<size_t, size_t> m_time_idx; // time, idx

            // Flush the block at the given index to disk, if it is dirty.
            void flushBlock(size_t idx);

            // Convert the row and column indices to an index.
            size_t toIdx(int32_t col, int32_t row);

            // Convert the index to a column index.
            int32_t toCol(size_t idx);

            // Convert the index to a row index.
            int32_t toRow(size_t idx);

            // Remove the oldest block from the cache and return a reference
            // to its memory. Flush first if dirty.
            T* freeOldest();

            // Remove a block from the cache and return a reference to its
            // memory for reuse. If there are more than the limit of blocks
            // keep removing until the cache is full-1.
            // Will flush dirty blocks.
            T* freeOne();

        public:
            BlockCache();

            // The number of columns in a single block.
            int32_t blockWidth();

            // The number of rows in a single block.
            int32_t blockHeight();

            // Convert the raster row and column into a row and
            // column index for a single block.
            size_t toBlockIdx(int32_t col, int32_t row);

            // Set the raster band for this cache to manage.
            void setRasterBand(GDALRasterBand *band);

            // Return true if the cache is managing a block that contains the
            // given column and row.
            bool hasBlock(int32_t col, int32_t row);

            // Return true if the given index is valid for this cache.
            bool hasBlock(size_t idx);

            // Set the number of blocks managed by the cache.
            void setSize(size_t size);

            // Return the number of blocks managed by the cache.
            size_t getSize();

            // Return a pointer to the block containing the given column and row.
            // If the block is to be written, the forWrite argument flags it as 
            // dirty.
            T* getBlock(int32_t col, int32_t row, bool forWrite);

            // Get the value from the cache.
            T get(int32_t col, int32_t row);

            // Set the value in the cache. Implies a dirty block.
            void set(int32_t col, int32_t row, T value);

            // Flush all blocks to disk.
            void flush();

            // Free all blocks.
            void close();

            ~BlockCache();

        };

        template <class T>
        class DLL_EXPORT Raster : public Grid<T> {
        private:
            int32_t m_cols, m_rows; // Raster cols/rows
            int32_t m_bandn; // The band number
            bool m_writable; // True if the raster is writable
            GDALDataset *m_ds; // GDAL dataset
            GDALRasterBand *m_band; // GDAL band
            GDALDataType m_type; // GDALDataType -- limits the possible template types.
            T m_nodata; // Nodata value.
            double m_trans[6]; // Raster transform
            bool m_inited; // True if the instance is initialized.
            std::string m_filename; // Raster filename
            BlockCache<T> m_cache; // Block cache.

            // Get the GDAL type for the given c++ type.
            GDALDataType getType(double v);

            GDALDataType getType(float v);

            GDALDataType getType(uint64_t v);

            GDALDataType getType(int64_t v);

            GDALDataType getType(uint32_t v);

            GDALDataType getType(int32_t v);

            GDALDataType getType(uint16_t v);

            GDALDataType getType(int16_t v);

            GDALDataType getType(uint8_t v);

            GDALDataType getType(int8_t v);

            // Get the GDAL type for the current raster.
            GDALDataType getType();

            // Return the default nodata for the raster's type.
            T getDefaultNodata();

        public:

            const static int32_t FLOAT64 = 1;
            const static int32_t FLOAT32 = 2;
            const static int32_t UINT32 = 3;
            const static int32_t UINT16 = 4;
            const static int32_t BYTE = 5;
            const static int32_t INT32 = 6;
            const static int32_t INT16 = 7;

            // Basic constructor.
            Raster();

            // Create a new raster for writing with a template of a different type.

            template <class U>
            Raster(const std::string &filename, int32_t band, const Raster<U> &tpl) {
                std::string proj;
                tpl.projection(proj);
                init(filename, band, tpl.minx(), tpl.miny(), tpl.maxx(), tpl.maxy(), tpl.resolutionX(),
                        tpl.resolutionY(), (T) tpl.nodata(), proj);
            }

            // Create a new raster for writing with a template of a different type.
            Raster(const std::string &filename, int32_t band, const Raster<T> &tpl);

            // Build a new raster with the given filename, bounds, resolution, nodata and projection.
            Raster(const std::string &filename, int32_t band, double minx, double miny, double maxx, double maxy,
                    double resolutionX, double resolutionY, double nodata, const std::string &proj);

            // Build a new raster with the given filename, bounds, resolution, nodata and projection.
            Raster(const std::string &filename, int32_t band, const Bounds &bounds, double resolutionX, double resolutionY, double nodata, int32_t crs);

            // Build a new raster with the given filename, bounds, resolution, nodata and SRID.
            Raster(const std::string &filename, int32_t band, double minx, double miny, double maxx, double maxy,
                    double resolutionX, double resolutionY, double nodata, int32_t crs);

            // Open the given raster and load the given band. Set the writable argument to true
            // to enable writing.
            Raster(const std::string &filename, int32_t band = 1, bool writable = false);

            // Initialize the raster using the given file (which may not exist) using
            // another raster as a template. The raster pixel types need not be the same.

            template <class U>
            void init(const std::string &filename, int32_t band, const Raster<U> &tpl) {
                g_debug("Raster init: " << filename << "; [tpl]");
                std::string proj;
                tpl.projection(proj);
                init(filename, band, tpl.minx(), tpl.miny(), tpl.maxx(), tpl.maxy(), tpl.resolutionX(),
                        tpl.resolutionY(), (double) tpl.nodata(), proj);
            }

            void init(const std::string &filename, int32_t band, const Bounds &bounds, double resolutionX, double resolutionY,
                    double nodata, const std::string &proj);

            // Initializes the raster with the given filename, bounds, resolution, nodata and projection.
            void init(const std::string &filename, int32_t band, double minx, double miny, double maxx, double maxy,
                    double resolutionX, double resolutionY, double nodata, const std::string &proj);

            // Initializes a Raster from the existing file.
            void init(const std::string &filename, int32_t band = 1, bool writable = false);

            // Attempts to return the datatype of the raster
            // with the given filename.
            static int32_t getType(const std::string &filename);

            // Set the number of blocks in the cache.
            void setCacheSize(size_t size);

            // Return the filename for this raster.
            std::string filename() const;

            // Return the number of bands in the raster.
            int32_t bandCount() const;

            // Converts a numerical (EPSG) crs code to a projection string.
            std::string epsg2ProjText(int32_t crs) const;

            // Set the band number. Invokes a flush and loads the new band.
            void setBand(int32_t band);

            // Returns the current band number.
            int32_t getBandNum();

            // Returns true if the raster is initialized.
            bool inited() const;

            void fill(const T value);

            void readBlock(int32_t col, int32_t row, Grid<T> &grd, int32_t dstCol = 0, int32_t dstRow = 0, int32_t xcols = 0, int32_t xrows = 0);

            void writeBlock(int32_t col, int32_t row, Grid<T> &grd, int32_t srcCol = 0, int32_t srcRow = 0, int32_t xcols = 0, int32_t xrows = 0);

            void writeBlock(Grid<T> &block);

            void readBlock(Grid<T> &block);

            // Get the x resolution.
            double resolutionX() const;

            // Get the y resolution.
            double resolutionY() const;

            bool positiveX() const;

            bool positiveY() const;

            // Write the projection data to the given string object.
            void projection(std::string &proj) const;

            // Return the GDAL datatype of the raster.
            GDALDataType type() const;

            // Return the raster's geographic bounds.
            geotools::util::Bounds bounds() const;

            // The minimum bounding x.
            double minx() const;

            // The maximum bounding x.
            double maxx() const;

            // The left x.
            double leftx() const;

            // The right x.
            double rightx() const;

            // The minimum bounding y.
            double miny() const;

            // The maximum bounding y.
            double maxy() const;

            // The top y.
            double topy() const;

            // The bottom y.
            double bottomy() const;

            // The width of the raster in map units.
            double width() const;

            // The height of the raster in map units.
            double height() const;

            T nodata() const;

            void nodata(const T nodata);

            int32_t cols() const;

            int32_t rows() const;

            // Returns the row for a given y-coordinate.
            int32_t toRow(double y) const;

            // Returns the column for a given x-coordinate.
            int32_t toCol(double x) const;

            // Returns the x-coordinate for a given column.
            double toX(int32_t col) const;

            // Returns the y-coordinate for a given row.
            double toY(int32_t row) const;

            // Returns the x-coordinate for the cell centroid of a given column.
            double toCentroidX(int32_t col) const;

            // Returns the y-coordinate for the cell centorid of a given row.
            double toCentroidY(int32_t row) const;

            size_t size() const;

            // Returns true if the pixel is nodata.
            bool isNoData(int32_t col, int32_t row);

            // Returns true if the pixel is nodata.
            bool isNoData(size_t idx);

            // Returns true if the pixel is nodata.
            bool isNoData(double x, double y);

            // Returns true if the pixel exists and is not nodata.
            bool isValid(int32_t c, int32_t r);

            // Returns true if the pixel exists and is not nodata.
            bool isValid(double x, double y);

            // Gets the pixel value or nodata if the pixel doesn't exist.
            T getOrNodata(double x, double y);

            // Gets the pixel value or nodata if the pixel doesn't exist.
            T getOrNodata(int32_t col, int32_t row);

            T *grid();

            bool hasGrid() const;

            T get(double x, double y);

            T get(int32_t col, int32_t row);

            T get(size_t idx);

            T operator[](size_t idx);

            void set(int32_t col, int32_t row, const T v);

            void set(size_t idx, const T v);

            void set(double x, double y, const T v);

            bool isSquare() const;

            bool has(int32_t col, int32_t row) const;

            bool has(double x, double y) const;

            bool has(size_t idx) const;

            // Flush the current block to the dataset.
            void flush();

            ~Raster();

        };

    } // raster

} // geotools


#endif
