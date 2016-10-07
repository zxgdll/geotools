#ifndef __SORTEDPOINTSTREAM_HPP__
#define __SORTEDPOINTSTREAM_HPP__

#include <map>
#include <string>
#include <fstream>
#include <memory>
#include <unordered_map>

#include "geos/geom/Geometry.h"

#include "util.hpp"

namespace geotools {
	namespace las {

		class LASPoint {
		public:
			static double scaleX, scaleY, scaleZ;
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

			void write(std::ostream &str) {
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


		};

		class SortedPointStream {
		private:
			unsigned int m_numBlocks;
			std::string m_file;
			std::map<unsigned int, std::unique_ptr<std::ifstream> > m_cacheFiles;
			std::unordered_map<unsigned int, unsigned int> m_cacheCounts;
			std::unordered_map<unsigned int, std::string> m_blockFilenames;
			geotools::util::Bounds m_fileBounds;
			geotools::util::Bounds m_currentBounds;
			unsigned int m_currentBlock;
			bool m_inited;
			unsigned int m_pointCount;

		public:

			/**
			 * Constructs a new point stream on the given file, with the given number of blocks.
			 * The mode indicates which dimensions will be stored.
			 * If mem is true, the dataset is stored in memory. If false, in files.
			 */
			SortedPointStream(const std::string &file, unsigned int numBlocks);

			~SortedPointStream();

			/**
			 * Initializes the point stream. Called implicitly by next()
			 * so should exit quickly.
			 */
			void init();

			geotools::util::Bounds fileBounds();

			geotools::util::Bounds currentBounds();

			unsigned int pointCount();

			/**
			 * Reads the next available point into the Point object.
			 * Updates the geometry with the extent of the read region (that is,
			 * the region within which there are no more points to read.)
			 * Returns false if there are no more points, true otherwise.
			 */
			bool next(LASPoint &pt, geos::geom::Geometry **geom);

		};

	}
}


#endif