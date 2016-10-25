#ifndef __SORTEDPOINTSTREAM_HPP__
#define __SORTEDPOINTSTREAM_HPP__

#include <map>
#include <string>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <mutex>
#include <thread>
#include <atomic>
#include <queue>

#include "liblas/liblas.hpp"

#include "geos/geom/Geometry.h"

#include "util.hpp"

using namespace geotools::util;

namespace geotools {
	namespace las {

		class LASPoint {
		public:
			static double scaleX, scaleY, scaleZ;
			const static unsigned long dataSize = 20;

			double x, y, z;
			uint16_t intensity;
			uint16_t ret, numRets;
			uint16_t scanDir;
			uint8_t cls;
			int8_t angle;

			LASPoint() {}
			LASPoint(const liblas::Point &pt) {
				update(pt);
			}

			static void setScale(double x, double y, double z) {
				scaleX = x;
				scaleY = y;
				scaleZ = z;
			}

			void update(const liblas::Point &pt) {
				x = pt.GetX();
				y = pt.GetY();
				z = pt.GetZ();
				intensity = pt.GetIntensity();
				ret = pt.GetReturnNumber();
				numRets = pt.GetNumberOfReturns();
				cls = pt.GetClassification().GetClass();
				angle = pt.GetScanAngleRank();
			}

			void read(std::istream &str) {
				int32_t xx, yy, zz;
				str.read((char *) &xx, sizeof(int32_t));
				str.read((char *) &yy, sizeof(int32_t));
				str.read((char *) &zz, sizeof(int32_t));
				str.read((char *) &intensity, sizeof(uint16_t));
				str.read((char *) &ret, sizeof(uint16_t));
				str.read((char *) &numRets, sizeof(uint16_t));
				str.read((char *) &cls, sizeof(uint8_t));
				str.read((char *) &angle, sizeof(int8_t));
				x = (double) (xx * scaleX);
				y = (double) (yy * scaleY);
				z = (double) (zz * scaleZ);
			}

			void write(std::ostream &str) const {
				int32_t xx = (int32_t) (x / scaleX);
				int32_t yy = (int32_t) (y / scaleY);
				int32_t zz = (int32_t) (z / scaleZ);
				str.write((char *) &xx, sizeof(int32_t));
				str.write((char *) &yy, sizeof(int32_t));
				str.write((char *) &zz, sizeof(int32_t));
				str.write((char *) &intensity, sizeof(uint16_t));
				str.write((char *) &ret, sizeof(uint16_t));
				str.write((char *) &numRets, sizeof(uint16_t));
				str.write((char *) &cls, sizeof(uint8_t));
				str.write((char *) &angle, sizeof(int8_t));
			}

			void write(std::FILE *str) const {
				int32_t xx = (int32_t) (x / scaleX);
				int32_t yy = (int32_t) (y / scaleY);
				int32_t zz = (int32_t) (z / scaleZ);
				std::fwrite((const void *) &xx, sizeof(int32_t), 1, str);
				std::fwrite((const void *) &yy, sizeof(int32_t), 1, str);
				std::fwrite((const void *) &zz, sizeof(int32_t), 1, str);
				std::fwrite((const void *) &intensity, sizeof(uint16_t), 1, str);
				std::fwrite((const void *) &ret, sizeof(uint16_t), 1, str);
				std::fwrite((const void *) &numRets, sizeof(uint16_t), 1, str);
				std::fwrite((const void *) &cls, sizeof(uint8_t), 1, str);
				std::fwrite((const void *) &angle, sizeof(int8_t), 1, str);
			}

			void read(std::FILE *str) {
				int32_t xx, yy, zz;
				if(!std::fread((void *) &xx, sizeof(int32_t), 1, str))
					g_runerr("Failed to read x.")
				if(!std::fread((void *) &yy, sizeof(int32_t), 1, str))
					g_runerr("Failed to read x.")
				if(!std::fread((void *) &zz, sizeof(int32_t), 1, str))
					g_runerr("Failed to read x.")
				if(!std::fread((void *) &intensity, sizeof(uint16_t), 1, str))
					g_runerr("Failed to read x.")
				if(!std::fread((void *) &ret, sizeof(uint16_t), 1, str))
					g_runerr("Failed to read x.")
				if(!std::fread((void *) &numRets, sizeof(uint16_t), 1, str))
					g_runerr("Failed to read x.")
				if(!std::fread((void *) &cls, sizeof(uint8_t), 1, str))
					g_runerr("Failed to read x.")
				if(!std::fread((void *) &angle, sizeof(int8_t), 1, str))
					g_runerr("Failed to read x.")
				x = (double) (xx * scaleX);
				y = (double) (yy * scaleY);
				z = (double) (zz * scaleZ);
			}

			bool last() const {
				return numRets > 0 && ret == numRets;
			}

			bool first() const {
				return numRets > 0 && ret == 1;
			}

			bool intermediate() const {
				return numRets > 2 && ret > 1 && ret < numRets;
			}
			
			bool ground() const {
				return cls == 2;
			}

			bool single() const {
				return numRets == 1;
			}
		};

		class PointFile {
		private:
			std::FILE *m_file;
			double m_blockSize;

			std::atomic_ulong m_pointCount;
			unsigned int m_nextRow;

			Bounds m_bounds;
			
			std::unordered_map<unsigned int, std::list<LASPoint> > m_cache;
			std::mutex m_mucache;

			std::unordered_map<unsigned int, unsigned int> m_row; // row jump
			std::map<unsigned int, unsigned int> m_firstRow; // row jump
			std::set<unsigned int> m_rows;
			std::string m_filename;
			
			std::queue<LASPoint> m_ptq;
			std::mutex m_muptq;

			std::atomic_bool m_running;
			std::list<std::thread> m_consumers;
			std::list<std::thread> m_writers;

			unsigned int getRow(const LASPoint &pt) {
				if(m_blockSize < 0) {
					return (unsigned int) ((pt.y - m_bounds.maxy()) / m_blockSize);
				} else {
					return (unsigned int) ((pt.y - m_bounds.miny()) / m_blockSize);
				}
			}

			unsigned long toRowOffset(unsigned int row) {
				return row * (sizeof(unsigned int) * 4 + LASPoint::dataSize * PointFile::rowLen);
			}

			unsigned long toPointOffset(unsigned int row, unsigned int pos) {
				return toRowOffset(row) + pos * LASPoint::dataSize;
			}

		public:

			const static unsigned int rowLen = 1024;

			// row: rowId | count | .. points .. 

			PointFile(const std::string &filename, const Bounds &bounds, double blockSize) {
				m_blockSize = blockSize;
				m_bounds.collapse();
				m_bounds.extend(bounds);
				m_filename = filename;
				m_file = nullptr;
			}

			bool running() {
				return m_running;
			}

			static void _consumer(PointFile *pf) {
				while(pf->running())
					pf->consumePoints();
			}

			static void _writer(PointFile *pf) {
				while(pf->running())
					pf->flushPoints();
			}

			void openWrite() {
				m_running = true;
				m_pointCount = 0;
				m_nextRow = 0;
				close();
				m_file = std::fopen(m_filename.c_str(), "wb");
				if(!m_file)
					g_runerr("Failed to open file: " << errno);
				for(int i = 0; i < 1; ++i) {
					m_consumers.push_back(std::thread(PointFile::_consumer, this));
					m_writers.push_back(std::thread(PointFile::_writer, this));
				}
			}

			void openRead() {
				m_nextRow = 0;
				close();
				m_file = std::fopen(m_filename.c_str(), "rb");
				if(!m_file)
					g_runerr("Failed to open file: " << errno);
			}

			void addPoint(const LASPoint &pt) {
				m_muptq.lock();
				m_ptq.push(pt);
				m_muptq.unlock();
			}

			void consumePoints() {
				if(m_ptq.size()) {
					m_muptq.lock();
					const LASPoint pt = m_ptq.front();
					m_ptq.pop();
					m_muptq.unlock();
					unsigned int row = getRow(pt);
					m_mucache.lock();
					m_cache[row].push_back(LASPoint(pt));
					m_rows.insert(row);
					m_mucache.unlock();
					++m_pointCount;
				}
			}

			void flushPoints() {
				for(const auto &it : m_cache) {
					if(it.second.size() >= PointFile::rowLen) {
						m_mucache.lock();
						unsigned int row = it.first;
						unsigned int nextRow = m_nextRow++;
						unsigned int lastRow = 0;
						unsigned int count = g_min(PointFile::rowLen, m_cache[row].size());
						if(m_row.find(row) != m_row.end())
							lastRow = m_row[row];
						auto it = m_cache[row].begin();
						auto it2 = it;
						std::advance(it2, count);
						std::vector<LASPoint> pts(it, it2);
						m_cache[row].erase(it, it2);
						flush(row, nextRow, lastRow, count, pts);
						m_row[row] = nextRow;
						m_mucache.unlock();
					}
				}
			}

			void flush(unsigned int row, unsigned int nextRow, unsigned int lastRow, unsigned int count, const std::vector<LASPoint> &pts) {

				unsigned int zero = 0;
				unsigned int tag = 1111;

				if(m_firstRow.find(row) == m_firstRow.end())
					m_firstRow[row] = nextRow;

				if(lastRow > 0) {
					//g_debug(" -- updating jump " << row << ", " << lastRow << ", " << nextRow);
					std::fseek(m_file, toRowOffset(lastRow) + 3 * sizeof(unsigned int), SEEK_SET);	
					std::fwrite((const void *) &nextRow, sizeof(unsigned int), 1, m_file);	
				}

				std::fseek(m_file, toRowOffset(nextRow), SEEK_SET);
				std::fwrite((const void *) &tag, sizeof(unsigned int), 1, m_file);
				std::fwrite((const void *) &row, sizeof(unsigned int), 1, m_file);
				std::fwrite((const void *) &count, sizeof(unsigned int), 1, m_file);
				std::fwrite((const void *) &zero, sizeof(unsigned int), 1, m_file);
				for(unsigned int i = 0; i < count; ++i) {
					const LASPoint &pt = pts[i];
					pt.write(m_file);
				}
		
			}

			bool nextRow(std::list<std::shared_ptr<LASPoint> > &pts) {
				while(m_rows.size()) {
					unsigned int row;
					bool cont = true;
					
					if(m_rows.size()) {						
						row = *(m_rows.begin());
						m_rows.erase(row);
						cont = false;
					}
						
					if(cont || m_firstRow.find(row) == m_firstRow.end())
						continue;

					unsigned int nextRow = m_firstRow[row];
					unsigned int count, curRow, tag;
					do {
						std::fseek(m_file, toRowOffset(nextRow), SEEK_SET);
						if(!std::fread((void *) &tag, sizeof(unsigned int), 1, m_file))
							g_runerr("Failed to read tag.");
						if(tag != 1111)
							g_runerr("Invalid row tag " << tag);
						if(!std::fread((void *) &curRow, sizeof(unsigned int), 1, m_file))
							g_runerr("Failed to read current row.");
						if(!std::fread((void *) &count, sizeof(unsigned int), 1, m_file))
							g_runerr("Failed to read count.");
						if(!std::fread((void *) &nextRow, sizeof(unsigned int), 1, m_file))
							g_runerr("Failed to read next row.");
						for(unsigned int i = 0; i < count; ++i) {
							std::shared_ptr<LASPoint> pt(new LASPoint());
							pt->read(m_file);
							pts.push_back(pt);
						}
						//g_debug(" -- jump " << row << ", " << nextRow);
					} while(nextRow > 0);
					return true;
				}
				return false;
			}

			void close() {
				if(m_consumers.size()) {
					m_running = false;
					for(std::thread &c : m_consumers)
						c.join();
					for(std::thread &w : m_writers)
						w.join();
					m_consumers.clear();
					m_writers.clear();
				}
				if(m_file) {
					std::fclose(m_file);
					m_file = nullptr;
				}
			}

			unsigned int pointCount() {
				return m_pointCount;
			}

			unsigned int rowCount() {
				return m_rows.size();
			}

			~PointFile() {
				close();
			}
		};



		class SortedPointStream {
		private:
			Bounds m_bounds;
			unsigned int m_pointCount;
			std::vector<unsigned int> m_rows;
			unsigned int m_row;
			PointFile *m_pf;

			void init(const std::list<std::string> &files, double blockSize, bool rebuild);

		public:

			/**
			 * Constructs a new point stream on the given files.
			 *
			 * The block size gives the size, in map units, of a block,
			 * blocks are square and stored as individual files.
			 */
			SortedPointStream(const std::list<std::string> &files, double blockSize, bool rebuild = true);

			~SortedPointStream();


			unsigned int pointCount() const;

			unsigned int rowCount() const;

			//bool contains(double x, double y) const;

			//bool contains(double x1, double y1, double x2, double y2) const;

			const Bounds& bounds() const;

			bool next(std::list<std::shared_ptr<LASPoint> > &pts);
		};

	}
}


#endif
