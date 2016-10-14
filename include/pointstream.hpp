#ifndef __POINTSTREAM_HPP__
#define __POINTSTREAM_HPP__

#include <string>
#include <fstream>

#include <liblas/liblas.hpp>

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

			LASPoint() : 
				ret(0), numRets(0) {}

			LASPoint(const liblas::Point &pt) : LASPoint() {
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

			bool last() {
				return numRets > 0 && ret == numRets;
			}

			bool first() {
				return numRets > 0 && ret == 1;
			}

			bool intermediate() {
				return numRets > 2 && ret > 1 && ret < numRets;
			}
			
			bool ground() {
				return cls == 2;
			}

			bool single() {
				return numRets == 1;
			}
		};

		class PointStream {
		private:
			unsigned int m_numBlocks;
			std::string m_file;
			geotools::util::Bounds m_fileBounds;
			unsigned int m_pointCount;
			liblas::Reader *m_lasReader;
			std::ifstream *m_instr;

			/**
			 * Initializes the point stream.
			 */
			void init(bool deepBounds);
		public:

			/**
			 * Constructs a new point stream on the given file.
			 * The mode indicates which dimensions will be stored.
			 */
			PointStream(const std::string &file, bool deepBounds = false);

			~PointStream();

			geotools::util::Bounds fileBounds();

			geotools::util::Bounds currentBounds();

			unsigned int pointCount();

			/**
			 * Reads the next available point into the LASPoint object.
			 */
			bool next(LASPoint &pt);

		};

	} // las

} // geotools


#endif