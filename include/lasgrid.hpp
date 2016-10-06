#ifndef __LASGRID_HPP__
#define __LASGRID_HPP__

#include <string>
#include <vector>
#include <set>
#include <map>
#include <memory>

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
#define TYPE_SKEW 12
#define TYPE_KURTOSIS 13
#define TYPE_RUGOSITY 14

#define LAS_EXT ".las"

#define ATT_HEIGHT 1
#define ATT_INTENSITY 2

namespace geotools {
	
	namespace las {

		namespace lasgrid_config {

			extern double defaultResolution;
			extern double defaultRadius;
			extern bool defaultSnapToGrid;
			extern unsigned char defaultType;
			extern unsigned char defaultAngleLimit;
			extern unsigned char defaultAttribute;
			extern std::set<unsigned char> defaultClasses;
			extern std::map<std::string, unsigned char> types;
			extern std::map<std::string, unsigned char> attributes;

		}

		/**
		 * Contains configuration values for running the lasgrid process.
		 */
		class LASGridConfig {
		public:
			std::string dstFile;
			std::list<std::string> lasFiles;
			std::set<unsigned char> classes;
			geotools::util::Bounds bounds;
			double radius;
			double resolution;
			unsigned short hsrid;
			unsigned short vsrid;
			unsigned char attribute;
			unsigned char type;
			unsigned char angleLimit;
			bool fill;
			bool snap;

			/**
			 * Interpret the attribute and  return the constant int value.
			 */
			unsigned char parseAtt(const std::string &attStr);

			/**
			 * Interpret the output type and return the constant int value.
			 */
			unsigned char parseType(const std::string &typeStr);
		};

		class LASGrid {
		private:
			std::shared_ptr<geotools::util::Callbacks> m_callbacks;
			std::shared_ptr<geotools::las::LASGridConfig> m_config;

			/**
			 * Check the configuration for validity. 
			 * Throw an exception if it's invalid or absent.
			 */
			void checkConfig(const LASGridConfig &config);

			/**
			 * Compute the working bounds and the selection of files
			 * to include in the working set.
			 */
			void computeWorkBounds(const std::list<std::string> &files, const Bounds &bounds, 
				std::set<std::string> &selectedFiles, Bounds &workBounds);

		public:
			/**
			 * Add a callbacks object to receive progress notifications.
			 */
			void setCallbacks(geotools::util::Callbacks *callbacks);

			/**
			 * Execute the gridding process.
			 */
			void lasgrid(const LASGridConfig &config);
		};

	} // las
	
} // geotools

#endif
