#ifndef __SORTEDPOINTSTREAM_HPP__
#define __SORTEDPOINTSTREAM_HPP__

#include <map>
#include <string>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
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
			~LASPoint();

			static void setScale(double x, double y, double z);

			void update(const liblas::Point &pt);

			void write(std::ostream &str) const;

			void read(std::istream &str);

			void write(std::FILE *str);

			void read(std::FILE *str);

			void write(void *str) const;

			void read(void *str);

			bool last() const;

			bool first() const;

			bool intermediate() const;
			
			bool ground() const;

			bool single() const;
		};


		class SortedPointStream {
		private:
			Bounds m_bounds;
			bool m_inited;
			bool m_rebuild;
			uint64_t m_pointCount;
			uint32_t m_rowCount;
			uint32_t m_row;
			uint32_t m_rowLen;
			uint32_t m_rowSize;
			uint64_t m_size;
			double m_blockSize;
			std::FILE *m_file;
			std::list<std::string> m_files;
			std::string m_filename;
			std::unordered_map<uint32_t, std::list<LASPoint*> > m_cache;
			std::unordered_map<uint32_t, uint32_t> m_jump;
			bool m_running;
			uint32_t m_nextJump;
			std::mutex m_fmtx;
			std::mutex m_cmtx;
			std::mutex m_wmtx;
			std::queue<std::string> m_fileq;
			std::unordered_set<uint32_t> m_flush;

		public:

			/**
			 * Constructs a new point stream on the given files.
			 *
			 * The block size gives the size, in map units, of a block,
			 * blocks are square and stored as individual files.
			 */
			SortedPointStream(const std::list<std::string> &files, double blockSize, bool rebuild = true);

			void init();

			~SortedPointStream();

			uint64_t pointCount() const;

			uint32_t rowCount() const;

			void produce();
			void consume();
			
			const Bounds& bounds() const;

			bool next(std::list<std::shared_ptr<LASPoint> > &pts);
		};

	}
}


#endif
