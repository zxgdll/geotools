#ifndef __LASGRID_HPP__
#define __LASGRID_HPP__

#include <string>
#include <vector>
#include <set>
#include <map>

#include "geotools.h"
#include "raster.hpp"
#include "util.hpp"

#define TYPE_MIN 1
#define TYPE_MAX 2
#define TYPE_MEAN 3
#define TYPE_DENSITY 4
#define TYPE_VARIANCE 7
#define TYPE_STDDEV 8
#define TYPE_PVARIANCE 12
#define TYPE_PSTDDEV 13
#define TYPE_COUNT 9
#define TYPE_QUANTILE 10
#define TYPE_MEDIAN 11

#define LAS_EXT ".las"

#define ATT_HEIGHT 1
#define ATT_INTENSITY 2

namespace geotools {
	
	namespace las {

		namespace lasgrid_config {

			extern double defaultResolution;
			extern double defaultRadius;
			extern bool defaultSnapToGrid;
			extern int defaultType;
			extern unsigned char defaultAngleLimit;
			extern std::set<int> defaultClasses;
			extern std::map<std::string, int> types;
			extern std::map<std::string, int> attributes;
			extern int defaultAttribute;

		}

		namespace lasgrid_util {

			int parseAtt(char *attStr);

			/**
			 * Interpret the output type and return the constant int value.
			 */
			int parseType(char *typeStr);

		}

		class LasGrid {
		private:
			void (*m_fileCallback)(float);
			void (*m_overallCallback)(float);
		public:
			void setFileCallback(void (*callback)(float));
			void setOverallCallback(void (*callback)(float));

			void lasgrid(std::string &dstFile, std::vector<std::string> &files, std::set<int> &classes,
				int crs, int attribute, int type, double radius,
				double resolution, geotools::util::Bounds &bounds, unsigned char angleLimit, bool fill, bool snap = true);
		};

	} // las
	
} // geotools

#endif
