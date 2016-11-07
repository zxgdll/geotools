#ifndef __LASGRID_HPP__
#define __LASGRID_HPP__

#include <string>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include <mutex>

#include "geotools.h"
#include "raster.hpp"
#include "util.hpp"
#include "laspoint.hpp"
#include "cellstats.hpp"
#include "finalizedpointstream.hpp"

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
#define TYPE_SKEW 16
#define TYPE_KURTOSIS 15
#define TYPE_RUGOSITY 14
#define TYPE_GAP_FRACTION 17

#define LAS_EXT ".las"

#define ATT_HEIGHT 1
#define ATT_INTENSITY 2

namespace geotools {
	
	namespace point {

		namespace pointstats_config {

			extern double defaultResolution;
			extern bool defaultSnapToGrid;
			extern bool defaultNormalize;
			extern unsigned char defaultType;
			extern unsigned char defaultAngleLimit;
			extern unsigned char defaultAttribute;
			extern unsigned char defaultGapFraction;
			extern unsigned int defaultQuantile;
			extern unsigned int defaultQuantiles;
			extern unsigned int defaultFilterQuantiles;
			extern unsigned int defaultFilterQuantileFrom;
			extern unsigned int defaultFilterQuantileTo;
			extern unsigned int defaultThreads;

			extern std::set<unsigned char> defaultClasses;
			extern std::map<std::string, unsigned char> types;
			extern std::map<std::string, unsigned char> attributes;
			extern std::map<std::string, unsigned char> gapFractionTypes;

		}

		/**
		 * Contains configuration values for running the lasgrid process.
		 */
		class PointStatsConfig {
		public:
			std::string dstFile;
			std::list<std::string> sourceFiles;
			std::set<unsigned char> classes;
			geotools::util::Bounds bounds;
			bool fill;
			bool snap;
			bool rebuild;
			bool normalize;
			double resolution;
			unsigned int threads;
			unsigned short hsrid;
			unsigned short vsrid;
			unsigned char attribute;
			unsigned char type;
			unsigned char angleLimit;
			unsigned char quantile;
			unsigned char quantiles;
			unsigned char gapFractionType;
			unsigned int quantileFilter;
			unsigned int quantileFilterFrom;
			unsigned int quantileFilterTo;

			/**
			 * Interpret the attribute and  return the constant int value.
			 */
			unsigned char parseAtt(const std::string &attStr);

			/**
			 * Interpret the output type and return the constant int value.
			 */
			unsigned char parseType(const std::string &typeStr);

			unsigned char parseGap(const std::string &typeStr);

			/**
			 * Returns true if the classes set contains the class.
			 */
			bool hasClass(unsigned char cls) const {
				return classes.find(cls) != classes.end();
			}

			bool hasClasses() const {
				return classes.size() > 0;
			}

			/**
			 * Returns true if the quantile filter does not pass all points.
			 */
			bool hasQuantileFilter() const {
				return quantileFilterFrom == 0 && quantileFilterTo == quantileFilter - 1;
			}

		};

		class PointStats {
		private:
			std::mutex m_cmtx;
			std::mutex m_qmtx;
			bool m_running;
			std::unordered_map<size_t, std::list<std::shared_ptr<geotools::las::LASPoint> > > m_cache;
			std::queue<size_t> m_idxq;
			std::unique_ptr<geotools::point::stats::CellStats> m_computer;
			std::unique_ptr<geotools::raster::MemRaster<float> > m_mem;
			std::unique_ptr<geotools::las::FinalizedPointStream> m_ps;

			/**
			 * Check the configuration for validity. 
			 * Throw an exception if it's invalid or absent.
			 */
			void checkConfig(const PointStatsConfig &config);

			geotools::point::stats::CellStats* getComputer(const PointStatsConfig &config);

			/**
			 * Compute the working bounds and the selection of files
			 * to include in the working set.
			 */
			void computeWorkBounds(const std::list<std::string> &files, const Bounds &bounds, 
				std::set<std::string> &selectedFiles, Bounds &workBounds, unsigned long *pointCount);

		public:

			/**
			 * Runs the compute loop. Inteded to be used by a thread.
			 */
			void runner();

			/**
			 * Runs the read loop. Inteded to be used by a thread.
			 */
			void reader();

			/**
			 * Execute the gridding process.
			 */
			void pointstats(const PointStatsConfig &config, const Callbacks *callbacks = nullptr);
		};

	} // las
	
} // geotools

#endif
