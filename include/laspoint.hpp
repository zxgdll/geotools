#ifndef __LASPOINT_HPP__
#define __LASPOINT_HPP__

#include <fstream>
#include <cstdio>

#include "liblas/liblas.hpp"

#include "util.hpp"

using namespace geotools::util;

namespace geotools {

	namespace las {

		class LASPoint {
		private:

			// Determine the scale factors used by all
			// points (coordinates are stored as ints in 
			// the LAS format.)
			//static double m_scaleX, m_scaleY, m_scaleZ;

			// Indicates the size taken up by a single
			// LASPoint in memory or on disc.
			const static uint64_t m_dataSize = 20;

		public:

			double x, y, z;
			uint16_t intensity;
			uint16_t returnNum, numReturns;
			uint16_t scanDirection;
			uint8_t cls;
			int8_t scanAngle;

			LASPoint();

			LASPoint(const liblas::Point &pt);

			~LASPoint();

			// Sets the scale values used by all points.
			static void setScale(double x, double y, double z);

			// Set fields from a liblas Point object.
			void update(const liblas::Point &pt);

			bool operator<(const LASPoint&) const;
			bool operator>(const LASPoint&) const;
			bool operator==(const LASPoint&) const;
			bool operator!=(const LASPoint&) const;
			bool operator<=(const LASPoint&) const;
			bool operator>=(const LASPoint&) const;

			void write(std::ostream &str) const;

			void read(std::istream &str);

			void write(std::FILE *str);

			void read(std::FILE *str);

			void write(void *str) const;

			void read(void *str);

			// Return true if this is a last return.
			bool last() const;

			// Return true if this is a first return.
			bool first() const;

			// Return true if this is an intermediate return.
			bool intermediate() const;
			
			// Return true if this is a ground point (class 2).
			bool ground() const;

			// Return true if this is a single return.
			bool single() const;

			static double scaleX();
			static double scaleY();
			static double scaleZ();
	
			static uint64_t dataSize();

		};
	}
}


#endif
