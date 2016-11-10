#ifndef __LASGRID_HPP__
#define __LASGRID_HPP__

#include <string>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include <mutex>
#include <condition_variable>

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
			extern uint8_t defaultType;
			extern uint8_t defaultAngleLimit;
			extern uint8_t defaultAttribute;
			extern uint8_t defaultGapFraction;
			extern uint32_t defaultQuantile;
			extern uint32_t defaultQuantiles;
			extern uint32_t defaultFilterQuantiles;
			extern uint32_t defaultFilterQuantileFrom;
			extern uint32_t defaultFilterQuantileTo;
			extern uint32_t defaultThreads;

			extern std::set<uint8_t> defaultClasses;
			extern std::map<std::string, uint8_t> types;
			extern std::map<std::string, uint8_t> attributes;
			extern std::map<std::string, uint8_t> gapFractionTypes;

		}

		/**
		 * Contains configuration values for running the lasgrid process.
		 */
		class PointStatsConfig {
		public:
			std::vector<std::string> dstFiles;
			std::vector<std::string> sourceFiles;
			std::vector<uint8_t> types;
			std::set<uint8_t> classes;
			geotools::util::Bounds bounds;
			bool fill;
			bool snap;
			bool rebuild;
			bool normalize;
			double resolution;
			uint32_t threads;
			uint16_t hsrid;
			uint16_t vsrid;
			uint8_t attribute;
			uint8_t angleLimit;
			uint8_t quantile;
			uint8_t quantiles;
			uint8_t gapFractionType;
			uint32_t quantileFilter;
			uint32_t quantileFilterFrom;
			uint32_t quantileFilterTo;

			/**
			 * Interpret the attribute and  return the constant int value.
			 */
			uint8_t parseAtt(const std::string &attStr);

			/**
			 * Interpret the output type and return the constant int value.
			 */
			std::vector<uint8_t> parseTypes(const std::vector<std::string> &typeStrs);

			uint8_t parseGap(const std::string &typeStr);

			/**
			 * Returns true if the classes set contains the class.
			 */
			bool hasClass(uint8_t cls) const {
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

		class blocking_queue {
		private:
			std::queue<size_t> m_q;
			std::condition_variable m_c;
			std::mutex m_m;
		public:
			bool pop(size_t *idx) {
				std::unique_lock<std::mutex> lk(m_m);
				m_c.wait(lk);
				if(m_q.empty())
					return false;
				*idx = m_q.front();
				m_q.pop();
				return true;
			}
			void push(size_t idx) {
				{
					std::lock_guard<std::mutex> lk(m_m);
					m_q.push(idx);
				}
				m_c.notify_one();
			}
			size_t size() {
				return m_q.size();
			}
			bool empty() {
				return m_q.empty();
			}
			void flush() {
				while(!empty())
					m_c.notify_one();
			}
			void finish() {
				m_c.notify_all();
			}
		};

		class PointStats {
		private:
			std::mutex m_cmtx;
			std::mutex m_qmtx;
			std::condition_variable m_cdn;

			bool m_running;
			std::unordered_map<size_t, std::list<std::shared_ptr<geotools::las::LASPoint> > > m_cache;
			std::queue<size_t> m_idxq;
			std::vector<std::unique_ptr<geotools::point::stats::CellStats> > m_computers;
			std::vector<std::unique_ptr<geotools::raster::MemRaster<float> > > m_mem;
			std::vector<std::unique_ptr<std::mutex> > m_mtx;
			std::unique_ptr<geotools::las::FinalizedPointStream> m_ps;
			blocking_queue m_bq;

			/**
			 * Check the configuration for validity. 
			 * Throw an exception if it's invalid or absent.
			 */
			void checkConfig(const PointStatsConfig &config);

			geotools::point::stats::CellStats* getComputer(const uint8_t &type, const PointStatsConfig &config);

			/**
			 * Compute the working bounds and the selection of files
			 * to include in the working set.
			 */
			void computeWorkBounds(const std::vector<std::string> &files, const Bounds &bounds, 
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
