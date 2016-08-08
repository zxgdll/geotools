#ifndef __LASGRID_HPP__
#define __LASGRID_HPP__

#include <string>
#include <vector>
#include <set>
#include <map>

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

		namespace config {

			double defaultResolution = 2.0;
			double defaultRadius = std::sqrt(_sq(defaultResolution / 2) * 2);
			bool defaultSnapToGrid = true;

			std::set<int> defaultClasses = {2};

			std::map<std::string, int> types = {
				{"Minimum", TYPE_MIN}, {"Maximum", TYPE_MAX}, {"Mean", TYPE_MEAN}, {"Density", TYPE_DENSITY},
				{"Sample Variance", TYPE_VARIANCE}, {"Sample Std. Dev.", TYPE_STDDEV}, {"Population Variance", TYPE_PVARIANCE},
				{"Population Std. Dev.", TYPE_PSTDDEV}, {"Count", TYPE_COUNT}, {"Quantile", TYPE_QUANTILE}, {"Median", TYPE_MEDIAN}
			};
			int defaultType = TYPE_MEAN;

			std::map<std::string, int> attributes {
				{"Height", ATT_HEIGHT}, {"Intensity", ATT_INTENSITY}
			};
			int defaultAttribute = ATT_HEIGHT;

		}

		namespace util {

			int parseAtt(char *attStr);

			/**
			 * Interpret the output type and return the constant int value.
			 */
			int parseType(char *typeStr);

		}

		void lasgrid(std::string &dstFile, std::vector<std::string> &files, std::set<int> &classes,
							int crs, int attribute, int type, double radius,
							double resolution, std::vector<double> &bounds, unsigned char angleLimit, bool fill);

	}
}

#endif
