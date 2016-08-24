#include <queue>

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
template <class U>
std::vector<int> Grid<T>::floodFill(int col, int row, FillOperator<T> &op, Grid<U> &other, U fill) {

	int minc = cols() + 1;
	int minr = rows() + 1;
	int maxc = -1;
	int maxr = -1;
	int area = 0;

	std::queue<std::unique_ptr<Cell> > q;
	q.push(std::unique_ptr<Cell>(new Cell(col, row)));

	std::vector<bool> visited(size(), false); // Tracks visited pixels.

	while(q.size()) {

		std::unique_ptr<Cell> cel = std::move(q.front());
		row = cel->row;
		col = cel->col;
		q.pop();

		size_t idx = (size_t) row * cols() + col;

		if(!visited[idx] && op.fill(get(col, row))) {

			minc = g_min(col, minc);
			maxc = g_max(col, maxc);
			minr = g_min(row, minr);
			maxr = g_max(row, maxr);
			++area;
			other.set(col, row, fill);
			visited[idx] = true;

			if(row > 0)
				q.push(std::unique_ptr<Cell>(new Cell(col, row - 1)));
			if(row < rows() - 1)
				q.push(std::unique_ptr<Cell>(new Cell(col, row + 1)));

			for(int c = col - 1; c >= 0; --c) {
				idx = (size_t) row * cols() + c;
				if(!visited[idx] && op.fill(get(c, row))) {
					minc = g_min(c, minc);
					++area;
					other.set(c, row, fill);
					visited[idx] = true;
					if(row > 0)
						q.push(std::unique_ptr<Cell>(new Cell(c, row - 1)));
					if(row < rows() - 1)
						q.push(std::unique_ptr<Cell>(new Cell(c, row + 1)));
				} else {
					break;
				}
			}
			for(int c = col + 1; c < cols(); ++c) {
				idx = (size_t) row * cols() + c;
				if(!visited[idx] && op.fill(get(c, row))) {
					maxc = g_max(c, maxc);
					++area;
					other.set(c, row, fill);
					visited[idx] = true;
					if(row > 0)
						q.push(std::unique_ptr<Cell>(new Cell(c, row - 1)));
					if(row < rows() - 1)
						q.push(std::unique_ptr<Cell>(new Cell(c, row + 1)));
				} else {
					break;
				}
			}
		}
	}
	return {minc, minr, maxc, maxr, area};
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
        g_trace("Smoothing grid...");
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
        int bufRows = g_max(1, g_min(rows(), (64 * 1024 * 1024) / sizeof(T) / cols()));
        g_trace(" - buffer rows: " << bufRows);

        double nd = nodata();
        size_t completed = 0;

        #pragma omp parallel for
        for(int row = 0; row < rows(); row += bufRows - size) {
                double weights[size * size];
                gaussianWeights(weights, size, sigma);
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
	m_nodata(0) {
}

template <class T>
MemRaster<T>::MemRaster(int cols, int rows) : MemRaster() {
	init(cols, rows);
}

template <class T>
template <class U>
MemRaster<T>::MemRaster(Grid<U> &tpl) : MemRaster() {
	init(tpl.cols(), tpl.rows());
}

template <class T>
MemRaster<T>::~MemRaster() {
	if(m_item_dealloc != nullptr) {
		for(size_t i = 0; i < (size_t) m_cols * m_rows; ++i)
			m_item_dealloc(m_grid[i]);
	}
	if(m_grid != nullptr)
		free(m_grid);
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
template <class U>
void MemRaster<T>::convert(MemRaster<U> &g) {
	g.init(cols(), rows());
	for(size_t i = 0; i < size(); ++i)
		g.set(i, (U) get(i));
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
void MemRaster<T>::init(int cols, int rows) {
	if(cols <= 0 || rows <= 0)
		g_argerr("Invalid row or column count.");
	if(cols != m_cols || rows != m_rows) {
		m_cols = cols;
		m_rows = rows;
		if(m_grid != nullptr)
			free(m_grid);
		m_grid = (T *) malloc(sizeof(T) * cols * rows);
	}
}

template <class T>
template <class U>
void MemRaster<T>::init(Grid<U> &tpl) {
	init(tpl.cols(), tpl.rows());
}

template <class T>
void MemRaster<T>::fill(const T value) {
	checkInit();
	for(size_t i = 0; i < size(); ++i)
		m_grid[i] = value;
}

template <class T>
T &MemRaster<T>::get(size_t idx) {
	checkInit();
	if(idx >= size())
		g_argerr("Index out of bounds.");
	return m_grid[idx];
}

template <class T>
T &MemRaster<T>::get(int col, int row) {
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
T &MemRaster<T>::operator[](size_t idx) {
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
void MemRaster<T>::readBlock(int col, int row, Grid<T> &block, int dstCol, int dstRow) {
	if(&block == this)
		g_argerr("Recursive call to readBlock.");
	if(col + block.cols() - dstRow > m_cols)
		g_argerr("Block is wider than the available space.");
	if(row + block.rows() - dstCol > m_rows)
		g_argerr("Block is taller than the available space.");
	if(block.hasGrid()) {
		for(int r = 0; r < block.rows() - dstRow; ++r) {
			std::memcpy(
				block.grid() + (r + dstRow) * block.cols() + dstCol,
				m_grid + (row + r) * m_cols + col,
				(block.cols() - dstCol) * sizeof(T)
			);
		}
	} else {
		for(int r = 0; r < block.rows() - dstRow; ++r) {
			for(int c = 0; c < block.cols() - dstCol; ++c)
				block.set(c + dstCol, r + dstRow, get(c + col, r + row));
		}
	}
}

template <class T>
void MemRaster<T>::writeBlock(int col, int row, Grid<T> &block, int srcCol, int srcRow) {
	if(&block == this)
		g_argerr("Recursive call to writeBlock.");
	if(col + block.cols() - srcCol > m_cols)
		g_argerr("Block is wider than the available space.");
	if(row + block.rows() - srcRow > m_rows)
		g_argerr("Block is taller than the available space.");
	if(block.hasGrid()) {
		for(int r = 0; r < block.rows(); ++r) {
			std::memcpy(
				m_grid + (row + r) * m_cols + col,
				block.grid() + (r + srcRow) * block.cols() + srcCol,
				(block.cols() - srcRow) * sizeof(T)
			);
		}
	} else {
		for(int r = 0; r < block.rows() - srcRow; ++r) {
			for(int c = 0; c < block.cols() - srcCol; ++c)
				set(c + col, r + row, block.get(c + srcCol, r + srcRow));
		}
	}
}

template <class T>
void MemRaster<T>::writeBlock(Grid<T> &block) {
	writeBlock(0, 0, block);
}

template <class T>
void MemRaster<T>::readBlock(Grid<T> &block) {
	readBlock(0, 0, block);
}


// Implementations for BlockCache
template <class T>
void BlockCache<T>::flushBlock(size_t idx) {
	if(hasBlock(idx) && m_band->GetDataset()->GetAccess() == GA_Update) {
		T *blk = m_blocks[idx];
		if(m_band->WriteBlock(toCol(idx) / m_bw, toRow(idx) / m_bh, blk) != CE_None)
			g_runerr("Failed to flush block.");
		m_dirty[idx] = false;
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
	auto it = m_time_idx.rbegin();
	size_t time = it->first;
	size_t idx = it->second;
	if(m_dirty[idx])
		flushBlock(idx);
	T *blk = m_blocks[idx];
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
	return m_blocks.find(idx) != m_blocks.end();
}

template <class T>
void BlockCache<T>::setSize(size_t size) {
	while(m_blocks.size() > size)
		freeOne();
	m_size = size;
}

template <class T>
size_t BlockCache<T>::getSize() {
	return m_size;
}

template <class T>
T* BlockCache<T>::getBlock(int col, int row, bool forWrite) {
	size_t idx = toIdx(col, row);
	if(!hasBlock(idx)) {
		T *blk = freeOne();
		if(!blk)
			blk = (T *) malloc(sizeof(T) * m_bw * m_bh);
		if(m_band->ReadBlock(col / m_bw, row / m_bh, blk) != CE_None)
			g_runerr("Failed to read block.");
		m_blocks[idx] = blk;
	}
	++m_time; // TODO: No provision for rollover
	m_time_idx.erase(m_idx_time[idx]);
	m_time_idx[m_time] = idx;
	m_idx_time[idx] = m_time;
	if(forWrite)
		m_dirty[idx] = true;
	return m_blocks[idx];
}

template <class T>
void BlockCache<T>::flush() {
	for(auto it = m_blocks.begin(); it != m_blocks.end(); ++it)
		flushBlock(it->first);
}

template <class T>
void BlockCache<T>::close() {
	flush();
	for(auto it = m_blocks.begin(); it != m_blocks.end(); ++it) {
		if(it->second)                  
			free(it->second);
	}
}

template <class T>
BlockCache<T>::~BlockCache() {
	close();
}


// Implementations for Raster
template <class T>
void Raster<T>::loadBlock(int col, int row, bool forWrite) {
	if(!m_inited)
		g_runerr("Not inited before attempted read.");
	if(!has(col, row))
		g_argerr("Row or column out of bounds:" << col << ", " << row);
	m_block = m_cache.getBlock(col, row, forWrite);
	if(!m_block)
		g_runerr("Failed to load block from cache.");
}

template <class T>
GDALDataType Raster<T>::getType(unsigned long v) {
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
	m_ds(nullptr), m_band(nullptr), m_block(nullptr),
	m_type(getType()) {
}

template <class T>
template <class U>
Raster<T>::Raster(const std::string &filename, int band, const Raster<U> &tpl) : Raster() {
	std::string proj;
	tpl.projection(proj);
	init(filename, band, tpl.minx(), tpl.miny(), tpl.maxx(), tpl.maxy(), tpl.resolutionX(),
		tpl.resolutionY(), (T) tpl.nodata(), proj);
}

//template <class T>
//Raster<T>::Raster(const std::string &filename, int band, const Raster<T> &tpl) : Raster() {
//	init(filename, band, tpl);
//}

template <class T>
Raster<T>::Raster(const std::string &filename, int band, double minx, double miny, double maxx, double maxy,
		double resolutionX, double resolutionY, double nodata, const std::string &proj) : Raster() {
	init(filename,band,  minx, miny, maxx, maxy, resolutionX, resolutionY, nodata, proj);
}

template <class T>
Raster<T>::Raster(const std::string &filename, int band, Bounds &bounds, double resolutionX, double resolutionY, double nodata, int crs) :
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
template <class U>
void Raster<T>::init(const std::string &filename, int band, const Raster<U> &tpl) {
	g_trace("Raster init: " << filename << "; [tpl]");
	std::string proj;
	tpl.projection(proj);
	init(filename, band, tpl.minx(), tpl.miny(), tpl.maxx(), tpl.maxy(), tpl.resolutionX(),
		tpl.resolutionY(), tpl.nodata(), proj);
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
	int width = (int) ((maxx - minx) / resolutionX) * (resolutionX < 0 ? -1 : 1);
	int height = (int) ((maxy - miny) / resolutionY) * (resolutionY < 0 ? -1 : 1);

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
void Raster<T>::setCacheSize(size_t size) {
	m_cache.setSize(size);
}

template <class T>
std::string Raster<T>::filename() const {
	return m_filename;
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
	MemRaster<T> grd(m_cache.blockWidth(), m_cache.blockHeight());
	grd.fill(value);
	for(int r = 0; r < rows() / grd.rows(); ++r) {
		for(int c = 0; c < cols() / grd.cols(); ++c) {
			if(m_band->WriteBlock(c, r, grd.grid()) != CE_None)
				g_runerr("Fill error.");
		}
	}
}

template <class T>
void Raster<T>::readBlock(int col, int row, Grid<T> &grd, int dstCol, int dstRow) {
	if(&grd == this)
		g_runerr("Recursive call to readBlock.");
	int cols = g_min(m_cols - col, grd.cols() - dstCol);
	int rows = g_min(m_rows - row, grd.rows() - dstRow);
	if(cols < 1 || rows < 1)
		g_argerr("Zero read size.");
	MemRaster<T> mr(cols, rows);
	if(m_band->RasterIO(GF_Read, col, row, cols, rows, mr.grid(), cols, rows, getType(), 0, 0) != CE_None)
		g_runerr("Failed to read from band.");
	grd.writeBlock(dstCol, dstRow, mr);
}

template <class T>
void Raster<T>::writeBlock(int col, int row, Grid<T> &grd, int srcCol, int srcRow) {
	if(&grd == this)
		g_runerr("Recursive call to writeBlock.");
	int cols = g_min(m_cols - col, grd.cols() - srcCol);
	int rows = g_min(m_rows - row, grd.rows() - srcRow);
	if(cols < 1 || rows < 1)
		g_argerr("Zero write size.");
	MemRaster<T> mr(cols, rows);
	grd.readBlock(srcCol, srcRow, mr);
	if(m_band->RasterIO(GF_Write, col, row, cols, rows, mr.grid(), cols, rows, getType(), 0, 0) != CE_None)
		g_runerr("Failed to write to band.");
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
void Raster<T>::projection(std::string &proj) const {
	proj.assign(m_ds->GetProjectionRef());
}

template <class T>
GDALDataType Raster<T>::type() const {
	return m_band->GetRasterDataType();
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
	return (int) ((y - m_trans[3]) / m_trans[5]);
}

template <class T>
int Raster<T>::toCol(double x) const {
	return (int) ((x - m_trans[0]) / m_trans[1]);
}

template <class T>
double Raster<T>::toX(int col) const {
	return (col * m_trans[1]) + m_trans[0];
}

template <class T>
double Raster<T>::toY(int row) const {
	return (row * m_trans[5]) + m_trans[3];
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
T &Raster<T>::get(double x, double y) {
	return get(toCol(x), toRow(y));
}

template <class T>
T &Raster<T>::get(int col, int row) {
	loadBlock(col, row, false);
	return m_block[m_cache.toBlockIdx(col, row)];
}

template <class T>
T &Raster<T>::get(size_t idx) {
	if(idx >= size())
		g_argerr("Index out of bounds.");
	return get(idx % m_cols, (int) idx / m_cols);
}

template <class T>
T &Raster<T>::operator[](size_t idx) {
	return get(idx);
}

template <class T>
void Raster<T>::set(int col, int row, T v) {
	if(!m_writable)
		g_runerr("This raster is not writable.");
	loadBlock(col, row, true);
	m_block[m_cache.toBlockIdx(col, row)] = v;
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


