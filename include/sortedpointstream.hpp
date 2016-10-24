#ifndef __SORTEDPOINTSTREAM_HPP__
#define __SORTEDPOINTSTREAM_HPP__

#include <map>
#include <string>
#include <fstream>
#include <memory>
#include <unordered_map>

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
				std::fwrite((const void *) &xx, sizeof(int32_t), 1, str); // 4
				std::fwrite((const void *) &yy, sizeof(int32_t), 1, str); // 4
				std::fwrite((const void *) &zz, sizeof(int32_t), 1, str); // 4
				std::fwrite((const void *) &intensity, sizeof(uint16_t), 1, str); // 2
				std::fwrite((const void *) &ret, sizeof(uint16_t), 1, str); // 2 
				std::fwrite((const void *) &numRets, sizeof(uint16_t), 1, str); // 2
				std::fwrite((const void *) &cls, sizeof(uint8_t), 1, str); // 1
				std::fwrite((const void *) &angle, sizeof(int8_t), 1, str); // 1
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
			unsigned int m_pointCount;
			Bounds m_bounds;
			std::unordered_map<unsigned int, std::list<LASPoint> > m_cache;
			std::unordered_map<unsigned int, unsigned int> m_row; // row jump
			std::map<unsigned int, unsigned int> m_firstRow; // row jump
			std::set<unsigned int> m_rows;
			unsigned int m_nextRow;
			std::string m_filename;

			unsigned int getRow(const LASPoint &pt) {
				if(m_blockSize < 0) {
					return (unsigned int) ((pt.y - m_bounds.maxy()) / m_blockSize);
				} else {
					return (unsigned int) ((pt.y - m_bounds.miny()) / m_blockSize);
				}
			}

			unsigned long toRowOffset(unsigned int row) {
				return row * (sizeof(unsigned int) * 2 + LASPoint::dataSize * PointFile::rowLen);
			}

			unsigned long toPointOffset(unsigned int row, unsigned int pos) {
				return toRowOffset(row) + pos * LASPoint::dataSize;
			}

		public:

			const static unsigned int rowLen = 5000;

			// row: rowId | count | .. points .. 

			PointFile(const std::string &filename, const Bounds &bounds, double blockSize) {
				m_blockSize = blockSize;
				m_bounds.collapse();
				m_bounds.extend(bounds);
				m_filename = filename;
				m_file = nullptr;
			}

			void openWrite() {
				m_pointCount = 0;
				m_nextRow = 0;
				close();
				m_file = std::fopen(m_filename.c_str(), "wb");
				if(!m_file)
					g_runerr("Failed to open file: " << errno);
			}

			void openRead() {
				m_nextRow = 0;
				close();
				m_file = std::fopen(m_filename.c_str(), "rb");
				if(!m_file)
					g_runerr("Failed to open file: " << errno);
			}

			void addPoint(const LASPoint &pt) {
				unsigned int row = getRow(pt);
				m_cache[row].push_back(LASPoint(pt));
				m_rows.insert(row);
				++m_pointCount;
				if(m_cache[row].size() >= PointFile::rowLen)
					flush(row);
			}

			void flush(unsigned int row) {
				unsigned int zero = 0;

				if(m_firstRow.find(row) == m_firstRow.end())
					m_firstRow[row] = m_nextRow;

				if(m_row.find(row) != m_row.end()) {
					std::fseek(m_file, toRowOffset(m_row[row]) + 2 * sizeof(unsigned int), SEEK_SET);	
					std::fwrite((const void *) &m_nextRow, sizeof(unsigned int), 1, m_file);	
				}
				unsigned int count = m_cache[row].size();
				std::fseek(m_file, toRowOffset(m_nextRow), SEEK_SET);
				std::fwrite((const void *) &row, sizeof(unsigned int), 1, m_file);
				std::fwrite((const void *) &count, sizeof(unsigned int), 1, m_file);
				std::fwrite((const void *) &zero, sizeof(unsigned int), 1, m_file);
				for(const LASPoint &pt : m_cache[row])
					pt.write(m_file);
				m_cache[row].clear();
				m_row[row] = m_nextRow;
				++m_nextRow;
			}

			void flushAll() {
				for(const auto &it : m_cache)
					flush(it.first);
			}

			bool nextRow(std::list<std::shared_ptr<LASPoint> > &pts) {
				if(!m_rows.size())
					return false;

				unsigned int row = *(m_rows.begin());
				m_rows.erase(row);

				unsigned int frow = m_firstRow[row];
				unsigned int count;
				unsigned int nextRow;
				do {
					std::fseek(m_file, toRowOffset(frow) + sizeof(unsigned int), SEEK_SET);
					std::fread((void *) &count, sizeof(unsigned int), 1, m_file);
					std::fread((void *) &nextRow, sizeof(unsigned int), 1, m_file);
					for(unsigned int i = 0; i < count; ++i) {
						std::shared_ptr<LASPoint> pt(new LASPoint());
						pt->read(m_file);
						pts.push_back(pt);
					}
				} while(nextRow > 0);
				return true;
			}

			void close() {
				if(m_file) {
					std::fclose(m_file);
					m_file = nullptr;
				}
			}

			unsigned int pointCount() {
				return m_pointCount;
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