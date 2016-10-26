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
			const static uint64_t dataSize = 20;

			double x, y, z;
			uint16_t intensity;
			uint16_t ret, numRets;
			uint16_t scanDir;
			uint8_t cls;
			int8_t angle;

			LASPoint();
			LASPoint(const liblas::Point &pt);

			static void setScale(double x, double y, double z);

			void update(const liblas::Point &pt);

			void write(std::ostream &str) const;

			void read(std::istream &str);

			void write(std::FILE *str) const;

			void read(std::FILE *str);

			void write(void *str) const;

			void read(void *str);

			bool last() const;

			bool first() const;

			bool intermediate() const;
			
			bool ground() const;

			bool single() const;
		};

		class PointFile {
		private:
			std::FILE *m_file;
			double m_blockSize;

			std::atomic_ulong m_pointCount;
			uint32_t m_nextRow;

			Bounds m_bounds;
			
			std::unordered_map<uint32_t, std::list<LASPoint> > m_cache;
			std::mutex m_mucache;

			std::unordered_map<uint32_t, uint32_t> m_row; // row jump
			std::map<uint32_t, uint32_t> m_firstRow; // row jump
			std::set<uint32_t> m_rows;
			std::string m_filename;
			
			std::queue<LASPoint> m_ptq;
			std::mutex m_muptq;

			std::atomic_bool m_running;
			std::list<std::thread> m_consumers;
			std::list<std::thread> m_writers;

			uint32_t getRow(const LASPoint &pt) {
				if(m_blockSize < 0) {
					return (uint32_t) ((pt.y - m_bounds.maxy()) / m_blockSize);
				} else {
					return (uint32_t) ((pt.y - m_bounds.miny()) / m_blockSize);
				}
			}

			uint64_t toRowOffset(uint32_t row) {
				return row * (sizeof(uint32_t) * 4 + LASPoint::dataSize * PointFile::rowLen);
			}

			uint64_t toPointOffset(uint32_t row, uint32_t pos) {
				return toRowOffset(row) + pos * LASPoint::dataSize;
			}

		public:

			const static uint32_t rowLen = 1024;

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
				for(uint32_t i = 0; i < 1; ++i) {
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
					uint32_t row = getRow(pt);
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
						uint32_t row = it.first;
						uint32_t nextRow = m_nextRow++;
						uint32_t lastRow = 0;
						uint32_t count = g_min(PointFile::rowLen, m_cache[row].size());
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

			void flush(uint32_t row, uint32_t nextRow, uint32_t lastRow, uint32_t count, const std::vector<LASPoint> &pts) {

				uint32_t zero = 0;
				uint32_t tag = 1111;

				if(m_firstRow.find(row) == m_firstRow.end())
					m_firstRow[row] = nextRow;

				if(lastRow > 0) {
					//g_debug(" -- updating jump " << row << ", " << lastRow << ", " << nextRow);
					std::fseek(m_file, toRowOffset(lastRow) + 3 * sizeof(uint32_t), SEEK_SET);	
					std::fwrite((const void *) &nextRow, sizeof(uint32_t), 1, m_file);	
				}

				std::fseek(m_file, toRowOffset(nextRow), SEEK_SET);
				std::fwrite((const void *) &tag, sizeof(uint32_t), 1, m_file);
				std::fwrite((const void *) &row, sizeof(uint32_t), 1, m_file);
				std::fwrite((const void *) &count, sizeof(uint32_t), 1, m_file);
				std::fwrite((const void *) &zero, sizeof(uint32_t), 1, m_file);
				for(uint32_t i = 0; i < count; ++i) {
					const LASPoint &pt = pts[i];
					pt.write(m_file);
				}
		
			}

			bool nextRow(std::list<std::shared_ptr<LASPoint> > &pts) {
				while(m_rows.size()) {
					uint32_t row;
					bool cont = true;
					
					if(m_rows.size()) {						
						row = *(m_rows.begin());
						m_rows.erase(row);
						cont = false;
					}
						
					if(cont || m_firstRow.find(row) == m_firstRow.end())
						continue;

					uint32_t nextRow = m_firstRow[row];
					uint32_t count, curRow, tag;
					do {
						std::fseek(m_file, toRowOffset(nextRow), SEEK_SET);
						if(!std::fread((void *) &tag, sizeof(uint32_t), 1, m_file))
							g_runerr("Failed to read tag.");
						if(tag != 1111)
							g_runerr("Invalid row tag " << tag);
						if(!std::fread((void *) &curRow, sizeof(uint32_t), 1, m_file))
							g_runerr("Failed to read current row.");
						if(!std::fread((void *) &count, sizeof(uint32_t), 1, m_file))
							g_runerr("Failed to read count.");
						if(!std::fread((void *) &nextRow, sizeof(uint32_t), 1, m_file))
							g_runerr("Failed to read next row.");
						for(uint32_t i = 0; i < count; ++i) {
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

			uint64_t pointCount() {
				return m_pointCount;
			}

			uint32_t rowCount() {
				return m_rows.size();
			}

			~PointFile() {
				close();
			}
		};



		class SortedPointStream {
		private:
			Bounds m_bounds;
			PointFile *m_pf;
			bool m_inited;
			bool m_rebuild;
			uint32_t m_pointCount;
			uint32_t m_row;
			double m_blockSize;
			std::list<std::string> m_files;
			std::vector<uint32_t> m_rows;

			void init();

		public:

			/**
			 * Constructs a new point stream on the given files.
			 *
			 * The block size gives the size, in map units, of a block,
			 * blocks are square and stored as individual files.
			 */
			SortedPointStream(const std::list<std::string> &files, double blockSize, bool rebuild = true);

			~SortedPointStream();

			uint64_t pointCount() const;

			uint32_t rowCount() const;

			const Bounds& bounds() const;

			bool next(std::list<std::shared_ptr<LASPoint> > &pts);
		};

	}
}


#endif
