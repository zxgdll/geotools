#ifndef __SORTEDPOINTSTREAM_HPP__
#define __SORTEDPOINTSTREAM_HPP__

#include <map>
#include <string>
#include <fstream>
#include <memory>

#include "geos/geom/Geometry.h"

#include "util.hpp"

namespace geotools {
	namespace las {

		class SortedPointStream {
		private:
			unsigned int m_numBlocks;
			unsigned char m_mode;
			bool m_mem;
			std::string m_file;
			std::map<unsigned int, std::unique_ptr<std::ifstream> > m_cacheFiles;
			std::map<unsigned int, unsigned int> m_cacheCounts;
			geotools::util::Bounds m_fileBounds;
			geotools::util::Bounds m_currentBounds;
			unsigned int m_currentBlock;
			bool m_inited;
			unsigned int m_pointCount;
			/**
			 * Initializes the point stream. Called implicitly by next()
			 * so should exit quickly.
			 */
			void init();

		public:
			static const unsigned char xy = 1;
			static const unsigned char z = 2;
			static const unsigned char xyz = 3;
			static const unsigned char intensity = 4;

			/**
			 * Constructs a new point stream on the given file, with the given number of blocks.
			 * The mode indicates which dimensions will be stored.
			 * If mem is true, the dataset is stored in memory. If false, in files.
			 */
			SortedPointStream(const std::string &file, unsigned int numBlocks, unsigned char mode, bool mem = false);

			~SortedPointStream();

			geotools::util::Bounds fileBounds();

			geotools::util::Bounds currentBounds();

			unsigned int pointCount();

			/**
			 * Reads the next available point into the Point object.
			 * Updates the geometry with the extent of the read region (that is,
			 * the region within which there are no more points to read.)
			 * Returns false if there are no more points, true otherwise.
			 */
			bool next(geotools::util::Point &pt, geos::geom::Geometry **geom);

		};

	}
}


#endif