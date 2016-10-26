#include <queue>
#include <string>
#include <fstream>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/filesystem.hpp>

#include "util.hpp"
#include "raster.hpp"

using namespace geotools::util;
using namespace geotools::raster;


// Implementations for Cell
Cell::Cell(int col, int row) :
	col(col), row(row) {
}


// Implementations for TargetOperator (for flood fill)
template <class T>
TargetOperator<T>::TargetOperator(T match) :
	m_match(match) {
}

template <class T>
bool TargetOperator<T>::fill(T value) const {
	return value == m_match;
}


// Implementations forthe Grid class
template <class T>
void Grid<T>::gaussianWeights(double *weights, int size, double sigma) {
	// If size is an even number, bump it up.
	if(size % 2 == 0) {
		++size;
		g_warn("Gaussian kernel size must be an odd number >=3. Bumping up to " << size);
	}
	for(int r = 0; r < size; ++r) {
        	for(int c = 0; c < size; ++c) {
			int x = size / 2 - c;
			int y = size / 2 - r;
			weights[r * size + c] = (1 / (2 * G_PI * sigma * sigma)) * pow(G_E, -((x * x + y * y) / (2.0 * sigma * sigma)));
		}
	}
}

template <class T>
void Grid<T>::computeStats() {
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
}

template <class T>
T Grid<T>::max() {
	if(!m_stats)
		computeStats();
	return m_max;
}

template <class T>
T Grid<T>::min() {
	if(!m_stats)
		computeStats();
	return m_min;
}

template <class T>
T Grid<T>::mean() {
	if(!m_stats)
		computeStats();
	return m_mean;
}

template <class T>
T Grid<T>::stddev() {
	if(!m_stats)
		computeStats();
	return m_stddev;
}

template <class T>
T Grid<T>::variance() {
	if(!m_stats)
		computeStats();
	return m_variance;
}

template <class T>
std::vector<int> Grid<T>::floodFill(int col, int row, T target, T fill) {
	TargetOperator<T> op(target);
	return floodFill(col, row, op, *this, fill);
}

template <class T>
std::vector<int> Grid<T>::floodFill(int col, int row, FillOperator<T> &op, T fill) {
	return floodFill(col, row, op, *this, fill);
}

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

template <class T>
void Grid<T>::smooth(Grid<T> &smoothed, double sigma, int size) {
        g_debug("Smoothing grid...");
        Util::status(0, 1, "Smoothing grid...");
        if(sigma <= 0)
                g_argerr("Sigma must be > 0.");
        if(size < 3)
                g_argerr("Kernel size must be 3 or larger.");
        if(size % 2 == 0) {
                g_warn("Kernel size must be odd. Rounding up.");
                size++;
        }

        // Guess at a good number of rows for each strip. Say 64MB each
        int bufRows = g_max((int) 1, g_min(rows(), (int) ((64 * 1024 * 1024) / sizeof(T) / cols())));
        g_debug(" - buffer rows: " << bufRows);

        double nd = nodata();
        size_t completed = 0;

        #pragma omp parallel for
        for(int row = 0; row < rows(); row += bufRows - size) {
                MemRaster<double> weights(size, size);
                gaussianWeights(weights.grid(), size, sigma);
                MemRaster<T> strip(cols(), g_min(bufRows, rows() - row));
                MemRaster<T> smooth(cols(), g_min(bufRows, rows() - row));
                MemRaster<T> buf(size, size);
                strip.nodata(nd);
                strip.fill(nd);
                smooth.nodata(nd);
                smooth.fill(nd);
                #pragma omp critical(a)
                {
                        // On the first loop, read from the first row and write to size/2 down
                        // On the other loops, read from row-size/2, and write to 0 down.
                        readBlock(0, row == 0 ? row : row - size / 2, strip, 0, row == 0 ? size / 2 : 0);
                }
                for(int r = 0; r < strip.rows() - size; ++r) {
                        for(int c = 0; c < strip.cols() - size; ++c) {
                                double v, t = 0.0;
                                bool foundNodata = false;
                                strip.readBlock(c, r, buf);
                                for(int gr = 0; !foundNodata && gr < size; ++gr) {
                                        for(int gc = 0; !foundNodata && gc < size; ++gc) {
                                                v = buf.get(gc, gr);
                                                if(v == nd) {
                                                        foundNodata = true;
                                                } else {
                                                        t += weights[gr * size + gc] * v;
                                                }
                                        }
                                }
                                if(!foundNodata)
                                        smooth.set(c + size / 2, r + size / 2, t);
                        }
                        #pragma omp atomic
                        ++completed;
                        Util::status(completed, rows(), "Smoothing grid...");
                }
                #pragma omp critical(b)
                {
                        // The blur buffer is always written size/2 down, so read from there.
                        smoothed.writeBlock(0, row, smooth, 0, size / 2);
                }
        }
        Util::status(1, 1, "Smoothing grid... Done.", true);

}


// Implementations for MemRaster
template <class T>
void MemRaster<T>::checkInit() const {
	if(m_grid == nullptr)
        	g_runerr("This instance has not been initialized.");
}

template <class T>
MemRaster<T>::MemRaster() :
	m_grid(nullptr),
	m_cols(-1), m_rows(-1),
	m_item_dealloc(nullptr),
	m_nodata(0),
	m_mmapped(false),
	m_size(0) {
}

template <class T>
MemRaster<T>::MemRaster(int cols, int rows, bool mapped) : MemRaster() {
	init(cols, rows, mapped);
}

template <class T>
MemRaster<T>::~MemRaster() {
	if(m_item_dealloc) {
		for(size_t i = 0; i < (size_t) m_cols * m_rows; ++i)
			m_item_dealloc(m_grid[i]);
	}
	freeMem();
}

template <class T>
void MemRaster<T>::setDeallocator(void (*item_dealloc)(T)) {
	m_item_dealloc = item_dealloc;
}

template <class T>
T *MemRaster<T>::grid() {
	return m_grid;
}

template <class T>
bool MemRaster<T>::hasGrid() const {
	return true;
}

template <class T>
int MemRaster<T>::rows() const {
	return m_rows;
}

template <class T>
int MemRaster<T>::cols() const {
	return m_cols;
}

template <class T>
size_t MemRaster<T>::size() const {
	return (size_t) m_rows * m_cols;
}

template <class T>
void MemRaster<T>::freeMem() {
	if(m_grid) {
		if(m_mmapped) {
			m_mappedFile.release();
		} else {
			free(m_grid);
		}
	}
}

template <class T>
void MemRaster<T>::init(int cols, int rows, bool mapped) {
	if(cols <= 0 || rows <= 0)
		g_argerr("Invalid row or column count.");
	if(cols != m_cols || rows != m_rows) {
		freeMem();
		m_cols = cols;
		m_rows = rows;
		m_mmapped = mapped;
		m_grid = nullptr;
		m_size = sizeof(T) * cols * rows;
		m_mappedFile.release();
		if(mapped) {
			const std::string filename = Util::tmpFile();
			m_mappedFile.reset(Util::mapFile(filename, m_size).release());
			m_grid = (T *) m_mappedFile->data();
		} else {
			m_grid = (T *) malloc(m_size);
		}
		if(!m_grid)
			g_runerr("Failed to allocate memory for MemRaster.");
	}
}

template <class T>
void MemRaster<T>::fill(const T value) {
	checkInit();
	for(size_t i = 0; i < size(); ++i)
		m_grid[i] = value;
}

template <class T>
T MemRaster<T>::get(size_t idx) {
	checkInit();
	if(idx >= size())
		g_argerr("Index out of bounds.");
	return m_grid[idx];
}

template <class T>
T MemRaster<T>::get(int col, int row) {
	size_t idx = (size_t) row * m_cols + col;
	return get(idx);
}

template <class T>
bool MemRaster<T>::isNoData(int col, int row) {
	return get(col, row) == m_nodata;
}

template <class T>
bool MemRaster<T>::isNoData(size_t idx) {
	return get(idx) == m_nodata;
}

template <class T>
void MemRaster<T>::set(int col, int row, const T value) {
	size_t idx = (size_t) row * m_cols + col;
	set(idx, value);
}

template <class T>
void MemRaster<T>::set(size_t idx, const T value) {
	checkInit();
	if(idx >= size())
		g_argerr("Index out of bounds.");
	m_grid[idx] = value;
}

template <class T>
bool MemRaster<T>::has(int col, int row) const {
	return col >= 0 && col < m_cols && row >= 0 && row < m_rows;
}

template <class T>
bool MemRaster<T>::has(size_t idx) const {
	return idx < (size_t) m_cols * m_rows;
}

template <class T>
T MemRaster<T>::operator[](size_t idx) {
	checkInit();
	if(idx >= size())
		g_argerr("Index out of bounds.");
	return m_grid[idx];
}

template <class T>
bool MemRaster<T>::isSquare() const {
	return cols() == rows();
}

template <class T>
void MemRaster<T>::toMatrix(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx) {
	for(int r = 1; r < rows(); ++r) {
		for(int c = 0; c < cols(); ++c)
			mtx(r, c) = get(c, r);
	}
}

template <class T>
void MemRaster<T>::fromMatrix(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx) {
	for(int r = 1; r < rows(); ++r) {
		for(int c = 0; c < cols(); ++c)
			set(c, r, mtx(r, c));
	}
}

template <class T>
T MemRaster<T>::nodata() const {
	return m_nodata;
}

template <class T>
void MemRaster<T>::nodata(T nodata) {
	m_nodata = nodata;
}

template <class T>
void MemRaster<T>::readBlock(int col, int row, Grid<T> &block, int dstCol, int dstRow, int xcols, int xrows) {
	if(&block == this)
		g_argerr("Recursive call to readBlock.");
	if(dstCol < 0 || dstRow < 0 || dstCol >= block.cols() || dstRow >= block.rows())
		g_argerr("Invalid destination column or row: row: " << dstRow << "; col: " << dstCol << "; block: " << block.rows() << "," << block.rows());
	int cols = g_min(m_cols - col, block.cols() - dstCol);
	int rows = g_min(m_rows - row, block.rows() - dstRow);
	if(block.hasGrid()) {
		// Copy rows of he source l
		for(int r = 0; r < rows; ++r) {
			std::memcpy(
				(block.grid() + (dstRow + r) * block.cols() + dstCol),
				(m_grid       + (row    + r) * m_cols       + col),
				cols * sizeof(T)
			);
		}
	} else {
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c)
				block.set(c + dstCol, r + dstRow, get(c + col, r + row));
		}
	}
}

template <class T>
void MemRaster<T>::writeBlock(int col, int row, Grid<T> &block, int srcCol, int srcRow, int xcols, int xrows) {
	if(&block == this)
		g_argerr("Recursive call to writeBlock.");
	if(srcCol < 0 || srcRow < 0 || srcCol >= block.cols() || srcRow >= block.rows())
		g_argerr("Invalid source column or row: row: " << srcRow << "; col: " << srcCol << "; block: " << block.rows() << "," << block.rows());
	int cols = g_min(m_cols - col, block.cols() - srcCol);
	int rows = g_min(m_rows - row, block.rows() - srcRow);
	if(xcols > 0) cols = g_min(xcols, cols);
	if(xrows > 0) rows = g_min(xrows, rows);
	if(block.hasGrid()) {
		for(int r = 0; r < rows; ++r) {
			std::memcpy(
				(m_grid       + (r + row)    * m_cols       + col),
				(block.grid() + (r + srcRow) * block.cols() + srcCol),
				cols * sizeof(T)
			);
		}
	} else {
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c)
				set(c + col, r + row, block.get(c + srcCol, r + srcRow));
		}
	}
}

template <class T>
void MemRaster<T>::writeBlock(Grid<T> &block) {
	writeBlock(0, 0, block, 0, 0, block.cols(), block.rows());
}

template <class T>
void MemRaster<T>::readBlock(Grid<T> &block) {
	readBlock(0, 0, block, 0, 0, block.cols(), block.rows());
}


// Implementations for BlockCache
template <class T>
void BlockCache<T>::flushBlock(size_t idx) {
	if(m_blocks.find(idx) != m_blocks.end() && m_band->GetDataset()->GetAccess() == GA_Update) {
		T *blk = m_blocks[idx];
		if(blk != nullptr) {
			#pragma omp critical(__gdal_io)
			{
				if(m_band->WriteBlock(toCol(idx) / m_bw, toRow(idx) / m_bh, blk) != CE_None)
					g_runerr("Failed to flush block.");
			}
			m_dirty[idx] = false;
		}
	}
}

template <class T>
size_t BlockCache<T>::toIdx(int col, int row) {
	return ((size_t) (col / m_bw) << 32) | (row / m_bh);
}

template <class T>
int BlockCache<T>::toCol(size_t idx) {
	return ((idx >> 32) & 0xffffffff) * m_bw;
}

template <class T>
int BlockCache<T>::toRow(size_t idx) {
	return (idx & 0xffffffff) * m_bh;
}

template <class T>
T* BlockCache<T>::freeOldest() {
	T *blk = nullptr;
	auto it = m_time_idx.rbegin();
	size_t time = it->first;
	size_t idx = it->second;
	if(m_dirty[idx])
		flushBlock(idx);
	blk = m_blocks[idx];
	m_blocks.erase(idx);
	m_idx_time.erase(idx);
	m_time_idx.erase(time);
	m_dirty.erase(idx);
	return blk;
}

template <class T>
T* BlockCache<T>::freeOne() {
	T *blk = nullptr;
	while(m_blocks.size() >= m_size) {
		if(blk)
			free(blk);
		blk = freeOldest();
	}
	return blk;
}

template <class T>
BlockCache<T>::BlockCache() :
	m_band(nullptr), 
	m_size(0),
	m_time(0),
	m_bw(0), m_bh(0) {
}

template <class T>
int BlockCache<T>::blockWidth() {
	return m_bw;
}

template <class T>
int BlockCache<T>::blockHeight() {
	return m_bh;
}

template <class T>
size_t BlockCache<T>::toBlockIdx(int col, int row) {
	return (row % m_bh) * m_bw + (col % m_bw);
}

template <class T>
void BlockCache<T>::setRasterBand(GDALRasterBand *band) {
	m_band = band;
	band->GetBlockSize(&m_bw, &m_bh);
}

template <class T>
bool BlockCache<T>::hasBlock(int col, int row) {
	return hasBlock(toIdx(col, row));
}

template <class T>
bool BlockCache<T>::hasBlock(size_t idx) {
	bool has = false;
	#pragma omp critical(__blockcache)
	{
		has = m_blocks.find(idx) != m_blocks.end();
	}
	return has;
}

template <class T>
void BlockCache<T>::setSize(size_t size) {
	#pragma omp critical(__blockcache)
	{
		while(m_blocks.size() > size)
			freeOne();
	}
	m_size = size;
}

template <class T>
size_t BlockCache<T>::getSize() {
	return m_size;
}

template <class T>
T* BlockCache<T>::getBlock(int col, int row, bool forWrite) {
	size_t idx = toIdx(col, row);
	T *blk = nullptr;
	if(m_blocks.find(idx) == m_blocks.end()) {
		//g_debug(" -- cache - new block");
		T *blk = freeOne();
		if(!blk)
			blk = (T *) malloc(sizeof(T) * m_bw * m_bh);
		if(!blk)
			g_runerr("Failed to allocate memory for raster block.");
		#pragma omp critical(__gdal_io)
		{
			if(m_band->ReadBlock(col / m_bw, row / m_bh, blk) != CE_None)
				g_runerr("Failed to read block.");
		}
		m_blocks[idx] = blk;
	}
	//g_debug(" -- cache - update time");
	++m_time; // TODO: No provision for rollover
	m_time_idx.erase(m_idx_time[idx]);
	m_time_idx[m_time] = idx;
	m_idx_time[idx] = m_time;
	if(forWrite)
		m_dirty[idx] = true;
	blk = m_blocks[idx];
	//g_debug(" -- cache - return block");
	return blk;
}

template <class T>
T BlockCache<T>::get(int col, int row) {
	T v;
	#pragma omp critical(__blockcache)
	{
		T *blk = getBlock(col, row, false);
		v = blk[toBlockIdx(col, row)];
	}
	return v;
}

template <class T>
void BlockCache<T>::set(int col, int row, T value) {
	#pragma omp critical(__blockcache)
	{	
		//g_debug(" -- cache set " << col << "," << row << "; " << value);
		T *blk = getBlock(col, row, true);
		blk[toBlockIdx(col, row)] = value;
	}
}

template <class T>
void BlockCache<T>::flush() {
	#pragma omp critical(__blockcache)
	{
		for(auto it = m_blocks.begin(); it != m_blocks.end(); ++it)
			flushBlock(it->first);
	}
}

template <class T>
void BlockCache<T>::close() {
	flush();
	#pragma omp critical(__blockcache)
	{
		for(auto it = m_blocks.begin(); it != m_blocks.end(); ++it) {
			if(it->second != nullptr) {
				free(it->second);
				it->second = nullptr;
			}
		}
	}
}

template <class T>
BlockCache<T>::~BlockCache() {
	close();
}


// Implementations for Raster
template <class T>
GDALDataType Raster<T>::getType(unsigned long v) {
	(void) v;
	g_runerr("Raster with 64 bit integral type requested. Not implemented.");
}

template <class T>
GDALDataType Raster<T>::getType(long v) {
	(void) v;
	g_runerr("Raster with 64 bit integral type requested. Not implemented.");
}

template <class T>
GDALDataType Raster<T>::getType(double v) {
	(void) v;
	return GDT_Float64;
}

template <class T>
GDALDataType Raster<T>::getType(float v) {
	(void) v;
	return GDT_Float32;
}

template <class T>
GDALDataType Raster<T>::getType(unsigned int v) {
	(void) v;
	return GDT_UInt32;
}

template <class T>
GDALDataType Raster<T>::getType(int v) {
	(void) v;
	return GDT_Int32;
}

template <class T>
GDALDataType Raster<T>::getType(unsigned short v) {
	(void) v;
	return GDT_UInt16;
}

template <class T>
GDALDataType Raster<T>::getType(short v) {
	(void) v;
	return GDT_Int16;
}

template <class T>
GDALDataType Raster<T>::getType(unsigned char v) {
	(void) v;
	return GDT_Byte;
}

template <class T>
GDALDataType Raster<T>::getType(char v) {
	(void) v;
	return GDT_Byte;
}

template <class T>
GDALDataType Raster<T>::getType() {
	return getType((T) 0);
}

template <class T>
T Raster<T>::getDefaultNodata() {
	switch(getType()) {
	case GDT_Float32:
	case GDT_Float64:
		return (T) -9999.0; // Hides the implicit overflow for unsigned types. (Compiler error.)
	default:
		return (T) 0;
	}
}

template <class T>
Raster<T>::Raster() :
	m_cols(-1), m_rows(-1),
	m_bandn(1),
	m_writable(false),
	m_ds(nullptr), m_band(nullptr),
	m_type(getType()) {
}

template <class T>
Raster<T>::Raster(const std::string &filename, int band, const Raster<T> &tpl) : Raster() {
	init(filename, band, tpl);
}

template <class T>
Raster<T>::Raster(const std::string &filename, int band, double minx, double miny, double maxx, double maxy,
		double resolutionX, double resolutionY, double nodata, const std::string &proj) : Raster() {
	init(filename,band,  minx, miny, maxx, maxy, resolutionX, resolutionY, nodata, proj);
}

template <class T>
Raster<T>::Raster(const std::string &filename, int band, const Bounds &bounds, double resolutionX, double resolutionY, double nodata, int crs) :
	Raster(filename, band, bounds.minx(), bounds.miny(), bounds.maxx(), bounds.maxy(),
	resolutionX, resolutionY, nodata, crs) {
}

template <class T>
Raster<T>::Raster(const std::string &filename, int band, double minx, double miny, double maxx, double maxy,
		double resolutionX, double resolutionY, double nodata, int crs) : Raster() {
	std::string proj = epsg2ProjText(crs);
	init(filename, band, minx, miny, maxx, maxy, resolutionX, resolutionY, nodata, proj);
}

template <class T>
Raster<T>::Raster(const std::string &filename, int band, bool writable) : Raster() {
	init(filename, band, writable);
}

template <class T>
void Raster<T>::init(const std::string &filename, int band, const Bounds &bounds, double resolutionX, double resolutionY,
	double nodata, const std::string &proj) {
	init(filename, band, bounds.minx(), bounds.miny(), bounds.maxx(), bounds.maxy(),
		resolutionX, resolutionY, nodata, proj);
}

template <class T>
void Raster<T>::init(const std::string &filename, int band, double minx, double miny, double maxx, double maxy,
		double resolutionX, double resolutionY, double nodata, const std::string &proj) {

	g_debug("Raster init: " << filename << ", " << minx << ", " << miny << ", " << maxx << ", " << maxy << ", " << resolutionX << ", " << resolutionY << ", " << nodata << ", " << proj);

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
	int width = (int) ((maxx - minx) / resolutionX) * (resolutionX < 0 ? -1 : 1) + 1;
	int height = (int) ((maxy - miny) / resolutionY) * (resolutionY < 0 ? -1 : 1) + 1;

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

template <class T>
void Raster<T>::init(const std::string &filename, int band, bool writable) {

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

template <class T>
int Raster<T>::getType(const std::string &filename) {
	GDALDataset *ds = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);
	int type = ds->GetRasterBand(1)->GetRasterDataType();
	GDALClose(ds);
	switch(type) {
	case GDT_Byte: return Raster<T>::BYTE;
	case GDT_UInt16: return Raster<T>::UINT16;
	case GDT_UInt32: return Raster<T>::UINT32;
	case GDT_Int16: return Raster<T>::INT16;
	case GDT_Int32: return Raster<T>::INT32;
	case GDT_Float32: return Raster<T>::FLOAT32;
	case GDT_Float64: return Raster<T>::FLOAT64;
	default:
		g_argerr("Unknown data type: " << type);
	}
}

template <class T>
void Raster<T>::setCacheSize(size_t size) {
	m_cache.setSize(size);
}

template <class T>
std::string Raster<T>::filename() const {
	return m_filename;
}

template <class T>
int Raster<T>::bandCount() const {
	return m_ds->GetRasterCount();
}

template <class T>
std::string Raster<T>::epsg2ProjText(int crs) const {
	OGRSpatialReference ref;
	char *wkt;
	ref.importFromEPSG(crs);
	ref.exportToWkt(&wkt);
	return std::string(wkt);
}

template <class T>
bool Raster<T>::inited() const {
	return m_inited;
}

template <class T>
void Raster<T>::fill(T value) {
	m_cache.flush();
	MemRaster<T> grd(m_cache.blockWidth(), m_cache.blockHeight());
	grd.fill(value);
	#pragma omp critical(__gdal_io)
	{
		for(int r = 0; r < rows() / grd.rows(); ++r) {
			for(int c = 0; c < cols() / grd.cols(); ++c) {
				if(m_band->WriteBlock(c, r, grd.grid()) != CE_None)
					g_runerr("Fill error.");
			}
		}
	}
}

template <class T>
void Raster<T>::readBlock(int col, int row, Grid<T> &grd, int dstCol, int dstRow, int xcols, int xrows) {
	if(&grd == this)
		g_runerr("Recursive call to readBlock.");
	m_cache.flush();
	int cols = g_min(m_cols - col, grd.cols() - dstCol);
	int rows = g_min(m_rows - row, grd.rows() - dstRow);
	if(xcols > 0) cols = g_min(xcols, cols);
	if(xrows > 0) rows = g_min(xrows, rows);
	if(cols < 1 || rows < 1)
		g_argerr("Zero read size.");
	if(grd.hasGrid()) {
		if(cols != grd.cols() || rows != grd.rows()) {
			MemRaster<T> g(cols, rows);
			#pragma omp critical(__gdal_io)
			{
				if(m_band->RasterIO(GF_Read, col, row, cols, rows, g.grid(), cols, rows, getType(), 0, 0) != CE_None)
					g_runerr("Failed to read from band. (1)");
				m_band->FlushCache();
			}
			grd.writeBlock(dstCol, dstRow, g);
		} else {
			#pragma omp critical(__gdal_io)
			{
				if(m_band->RasterIO(GF_Read, col, row, cols, rows, grd.grid(), cols, rows, getType(), 0, 0) != CE_None)
					g_runerr("Failed to read from band. (2)");
				m_band->FlushCache();
			}
		}
	} else {
		MemRaster<T> mr(cols, rows);
		#pragma omp critical(__gdal_io)
		{
			if(m_band->RasterIO(GF_Read, col, row, cols, rows, mr.grid(), cols, rows, getType(), 0, 0) != CE_None)
				g_runerr("Failed to read from band. (3)");
			m_band->FlushCache();
		}
		grd.writeBlock(dstCol, dstRow, mr);
	}
}

template <class T>
void Raster<T>::writeBlock(int col, int row, Grid<T> &grd, int srcCol, int srcRow, int xcols, int xrows) {
	if(&grd == this)
		g_runerr("Recursive call to writeBlock.");
	m_cache.flush();
	int cols = g_min(m_cols - col, grd.cols() - srcCol);
	int rows = g_min(m_rows - row, grd.rows() - srcRow);
	if(xcols > 0) cols = g_min(xcols, cols);
	if(xrows > 0) rows = g_min(xrows, rows);
	if(cols < 1 || rows < 1)
		g_argerr("Zero write size.");
	if(grd.hasGrid()) {
		if(cols != grd.cols() || rows != grd.rows()) {
			MemRaster<T> g(cols, rows);
			grd.readBlock(srcCol, srcRow, g);
			#pragma omp critical(__gdal_io)
			{
				if(m_band->RasterIO(GF_Write, col, row, cols, rows, g.grid(), cols, rows, getType(), 0, 0) != CE_None)
					g_runerr("Failed to write to band. (1)");
				m_band->FlushCache();
			}
		} else {
			#pragma omp critical(__gdal_io)
			{
				if(m_band->RasterIO(GF_Write, col, row, cols, rows, grd.grid(), cols, rows, getType(), 0, 0) != CE_None)
					g_runerr("Failed to write to band. (2)");
				m_band->FlushCache();
			}
		}
	} else {
		MemRaster<T> mr(cols, rows);
		grd.readBlock(srcCol, srcRow, mr);
		#pragma omp critical(__gdal_io)
		{
			if(m_band->RasterIO(GF_Write, col, row, cols, rows, mr.grid(), cols, rows, getType(), 0, 0) != CE_None)
				g_runerr("Failed to write to band. (3)");
			m_band->FlushCache();
		}
	}
}

template <class T>
void Raster<T>::writeBlock(Grid<T> &block) {
	writeBlock(0, 0, block);
}

template <class T>
void Raster<T>::readBlock(Grid<T> &block) {
	readBlock(0, 0, block);
}

template <class T>
double Raster<T>::resolutionX() const {
	return m_trans[1];
}

template <class T>
double Raster<T>::resolutionY() const {
	return m_trans[5];
}

template <class T>
bool Raster<T>::positiveX() const {
	return resolutionX() > 0;
}

template <class T>
bool Raster<T>::positiveY() const {
	return resolutionY() > 0;
}

template <class T>
void Raster<T>::setBand(int band) {
	if(band == m_bandn)
		return;
	flush();
	m_band = m_ds->GetRasterBand(band);
	m_bandn = band;
}
template <class T>
int Raster<T>::getBandNum() {
	return m_bandn;
}

template <class T>
void Raster<T>::projection(std::string &proj) const {
	proj.assign(m_ds->GetProjectionRef());
}

template <class T>
GDALDataType Raster<T>::type() const {
	return m_band->GetRasterDataType();
}

template <class T>
Bounds Raster<T>::bounds() const {
	return Bounds(minx(), miny(), maxx(), maxy());
}

template <class T>
double Raster<T>::minx() const {
	return toX(0);
}

template <class T>
double Raster<T>::maxx() const {
	return toX(cols());
}

template <class T>
double Raster<T>::miny() const {
	return toY(rows());
}

template <class T>
double Raster<T>::maxy() const {
	return toY(0);
}

template <class T>
double Raster<T>::leftx() const {
	return resolutionX() > 0 ? minx() : maxx();
}

template <class T>
double Raster<T>::rightx() const {
	return resolutionX() < 0 ? minx() : maxx();
}

template <class T>
double Raster<T>::topy() const {
	return resolutionY() > 0 ? miny() : maxy();
}

template <class T>
double Raster<T>::bottomy() const {
	return resolutionY() < 0 ? miny() : maxy();
}

template <class T>
double Raster<T>::width() const {
	return maxx() - minx();
}

template <class T>
double Raster<T>::height() const {
	return maxy() - miny();
}

template <class T>
T Raster<T>::nodata() const {
	return m_nodata;
}

template <class T>
void Raster<T>::nodata(T nodata) {
	m_band->SetNoDataValue(nodata);
	m_nodata = nodata;
}

template <class T>
int Raster<T>::cols() const {
	return m_cols;
}

template <class T>
int Raster<T>::rows() const {
	return m_rows;
}

template <class T>
int Raster<T>::toRow(double y) const {
	return (y - m_trans[3]) / m_trans[5];
}

template <class T>
int Raster<T>::toCol(double x) const {
	return (x - m_trans[0]) / m_trans[1];
}

template <class T>
double Raster<T>::toX(int col) const {
	return m_trans[0] + col * m_trans[1] + m_trans[1] / 2.0;
}

template <class T>
double Raster<T>::toY(int row) const {
	return m_trans[3] + row * m_trans[5] + m_trans[5] / 2.0;
}

template <class T>
size_t Raster<T>::size() const {
	return (size_t) (m_cols * m_rows);
}

template <class T>
bool Raster<T>::isNoData(int col, int row) {
	return get(col, row) == m_nodata;
}

template <class T>
bool Raster<T>::isNoData(size_t idx) {
	return get(idx) == m_nodata;
}

template <class T>
bool Raster<T>::isNoData(double x, double y) {
	return isNoData(toCol(x), toRow(y));
}

template <class T>
bool Raster<T>::isValid(int c, int r) {
	return getOrNodata(c, r) != m_nodata;
}

template <class T>
bool Raster<T>::isValid(double x, double y) {
	return getOrNodata(x, y) != m_nodata;
}

template <class T>
T Raster<T>::getOrNodata(double x, double y) {
	if(!has(x, y)) {
		return m_nodata;
	} else {
		return get(x, y);
	}
}

template <class T>
T Raster<T>::getOrNodata(int col, int row) {
	if(!has(col, row)) {
		return m_nodata;
	} else {
		return get(col, row);
	}
}

template <class T>
T *Raster<T>::grid() {
	g_implerr("grid() Not implemented in Raster.");
}

template <class T>
bool Raster<T>::hasGrid() const {
	return false;
}

template <class T>
T Raster<T>::get(double x, double y) {
	return get(toCol(x), toRow(y));
}

template <class T>
T Raster<T>::get(int col, int row) {
	return m_cache.get(col, row);
}

template <class T>
T Raster<T>::get(size_t idx) {
	if(idx >= size())
		g_argerr("Index out of bounds.");
	return get(idx % m_cols, (int) idx / m_cols);
}

template <class T>
T Raster<T>::operator[](size_t idx) {
	return get(idx);
}

template <class T>
void Raster<T>::set(int col, int row, T v) {
	if(!m_writable)
		g_runerr("This raster is not writable.");
	m_cache.set(col, row, v);
}

template <class T>
void Raster<T>::set(size_t idx, T v) {
	if(idx >= size())
		g_argerr("Index out of bounds.");
	set(idx % m_cols, (int) idx / m_cols, v);
}

template <class T>
void Raster<T>::set(double x, double y, T v) {
	set(toCol(x), toRow(y), v);
}

template <class T>
bool Raster<T>::isSquare() const {
	return cols() == rows();
}

template <class T>
bool Raster<T>::has(int col, int row) const {
	return col >= 0 && col < m_cols && row >= 0 && row < m_rows;
}

template <class T>
bool Raster<T>::has(double x, double y) const {
	return has(toCol(x), toRow(y));
}

template <class T>
bool Raster<T>::has(size_t idx) const {
	return idx < (size_t) (m_cols * m_rows);
}

template <class T>
void Raster<T>::flush() {
	if(m_writable)
		m_cache.flush();
}

template <class T>
Raster<T>::~Raster() {
	m_cache.close();
	if(m_ds) // Probably not necessary.
		GDALClose(m_ds);
}


template class geotools::raster::Grid<float>;
template class geotools::raster::Grid<double>;
template class geotools::raster::Grid<unsigned int>;
template class geotools::raster::Grid<unsigned short>;
template class geotools::raster::Grid<unsigned char>;
template class geotools::raster::Grid<unsigned long>;
template class geotools::raster::Grid<int>;
template class geotools::raster::Grid<short>;
template class geotools::raster::Grid<char>;
template class geotools::raster::Grid<long>;

template class geotools::raster::Raster<float>;
template class geotools::raster::Raster<double>;
template class geotools::raster::Raster<unsigned int>;
template class geotools::raster::Raster<unsigned short>;
template class geotools::raster::Raster<unsigned char>;
template class geotools::raster::Raster<unsigned long>;
template class geotools::raster::Raster<int>;
template class geotools::raster::Raster<short>;
template class geotools::raster::Raster<char>;
template class geotools::raster::Raster<long>;

template class geotools::raster::MemRaster<float>;
template class geotools::raster::MemRaster<double>;
template class geotools::raster::MemRaster<unsigned int>;
template class geotools::raster::MemRaster<unsigned short>;
template class geotools::raster::MemRaster<unsigned char>;
template class geotools::raster::MemRaster<unsigned long>;
template class geotools::raster::MemRaster<int>;
template class geotools::raster::MemRaster<short>;
template class geotools::raster::MemRaster<char>;
template class geotools::raster::MemRaster<long>;
template class geotools::raster::MemRaster<std::vector<double>*>;


