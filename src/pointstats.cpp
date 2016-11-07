/*
 * Grids a point cloud represented by one or more LAS files.
 * Can produce grids of from intensity and elevation, using
 * minimum, maximum, mean, std dev, density, variance and count.
 *
 * Authored by: Rob Skelly rob@dijital.ca
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <list>
#include <climits>
#include <memory>
#include <cstring>
#include <cstdio>
#include <math.h>
#include <exception>
#include <unordered_set>
#include <cmath>
#include <thread>
#include <mutex>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "pointstats.hpp"
#include "lasutil.hpp"
#include "laspoint.hpp"
#include "finalizedpointstream.hpp"
#include "cellstats.hpp"

using namespace geotools::util;
using namespace geotools::raster;
using namespace geotools::point;
using namespace geotools::las;
using namespace geotools::point::stats;

namespace geotools {

	namespace point {

		namespace pointstats_config {
			
			double defaultResolution = 10.0; 
			bool defaultSnapToGrid = true;
			bool defaultNormalize = false;
			uint8_t defaultType = TYPE_MEAN;
			uint8_t defaultAttribute = ATT_HEIGHT;
			uint8_t defaultAngleLimit = 180;
			uint8_t defaultGapFraction = GAP_BLB;
			uint32_t defaultQuantile = 49;
			uint32_t defaultQuantiles = 100;
			uint32_t defaultFilterQuantiles = 100;
			uint32_t defaultFilterQuantileFrom = 0;
			uint32_t defaultFilterQuantileTo = 99;
			uint32_t defaultThreads = 1;

			std::set<uint8_t> defaultClasses = {2}; 

			std::map<std::string, uint8_t> types = {
				{"Minimum", TYPE_MIN}, {"Maximum", TYPE_MAX}, {"Mean", TYPE_MEAN}, {"Density", TYPE_DENSITY},
				{"Sample Variance", TYPE_VARIANCE}, {"Sample Std. Dev.", TYPE_STDDEV}, {"Population Variance", TYPE_PVARIANCE},
				{"Population Std. Dev.", TYPE_PSTDDEV}, {"Count", TYPE_COUNT}, {"Quantile", TYPE_QUANTILE}, 
				{"Median", TYPE_MEDIAN}, {"Rugosity", TYPE_RUGOSITY}, {"Kurtosis", TYPE_KURTOSIS}, {"Skewness", TYPE_SKEW},
				{"Gap Fraction", TYPE_GAP_FRACTION}
			};
			std::map<std::string, uint8_t> attributes = {
				{"Height", ATT_HEIGHT}, {"Intensity", ATT_INTENSITY}
			};
			std::map<std::string, uint8_t> gapFractionTypes = {
				{"IR", GAP_IR}, {"BLa", GAP_BLA}, {"BLb", GAP_BLB}, {"RR", GAP_RR}, {"FR", GAP_FR}
			};

		} // config
		
		uint8_t PointStatsConfig::parseAtt(const std::string &attStr) {
			if("intensity" == attStr) {
				return ATT_INTENSITY;
			} else if("height" == attStr) {
				return ATT_HEIGHT;
			} 
			return 0;
		}

		uint8_t PointStatsConfig::parseType(const std::string &typeStr) {
			if("min" == typeStr) {
				return TYPE_MIN;
			} else if("max" == typeStr) {
				return TYPE_MAX;
			} else if("mean" == typeStr) {
				return TYPE_MEAN;
			} else if("density" == typeStr) {
				return TYPE_DENSITY;
			} else if("variance" == typeStr) {
				return TYPE_VARIANCE;
			} else if("stddev" == typeStr) {
				return TYPE_STDDEV;
			} else if("pvariance" == typeStr) {
				return TYPE_PVARIANCE;
			} else if("pstddev" == typeStr) {
				return TYPE_PSTDDEV;
			} else if("count" == typeStr) {
				return TYPE_COUNT;
			} else if("median" == typeStr) {
				return TYPE_MEDIAN;
			} else if("skew" == typeStr) {
				return TYPE_SKEW;
			} else if("rugosity" == typeStr) {
				return TYPE_RUGOSITY;
			} else if("kurtosis" == typeStr) {
				return TYPE_KURTOSIS;
			} else if("gap" == typeStr) {
				return TYPE_GAP_FRACTION;
			}
			return 0;
		}

		uint8_t PointStatsConfig::parseGap(const std::string &gapStr) {
			if("bla" == gapStr) {
				return GAP_BLA;
			} else if("blb" == gapStr) {
				return GAP_BLB;
			} else if("fr" == gapStr) {
				return GAP_FR;
			} else if("rr" == gapStr) {
				return GAP_RR;
			} else if("ir" == gapStr) {
				return GAP_IR;
			}
			return 0;
		}

		void PointStats::checkConfig(const PointStatsConfig &config) {
			if(config.resolution <= 0.0)
				g_argerr("Resolution must be > 0: " << config.resolution);
			if(config.sourceFiles.size() == 0)
				g_argerr("At least one input file is required.");
			if(config.dstFile.empty()) 
				g_argerr("An output file is required.");
			if(config.attribute == 0)
				g_argerr("An attribute is required.");
			if(config.type == 0)
				g_argerr("A valid type is required.");
			if(config.classes.size() == 0)
				g_warn("No classes given. Matching all classes.");
			if(config.angleLimit <= 0)
				g_argerr("Angle limit must be greater than zero.");

			g_debug("Resolution: " << config.resolution);
			g_debug("Files: " << config.sourceFiles.size());
			g_debug("Destination: " << config.dstFile);
			g_debug("Attribute: " << config.attribute);
			g_debug("Type: " << config.type);
			g_debug("Classes: " << config.classes.size());
			g_debug("Angle Limit: " << config.angleLimit);
		}

		CellStats* PointStats::getComputer(const PointStatsConfig &config) {
			using namespace geotools::point::stats;
			switch(config.type) {
			case TYPE_MEAN: 		return new CellMean();
			case TYPE_MEDIAN: 		return new CellMedian();
			case TYPE_COUNT: 		return new CellCount();
			case TYPE_STDDEV: 		return new CellSampleStdDev();
			case TYPE_VARIANCE: 	return new CellSampleVariance();
			case TYPE_PSTDDEV: 		return new CellPopulationStdDev();
			case TYPE_PVARIANCE:	return new CellPopulationVariance();
			case TYPE_DENSITY: 		return new CellDensity(g_sq(config.resolution));
			case TYPE_RUGOSITY: 	return new CellRugosity(config.resolution * config.resolution, 20.0);
			case TYPE_MAX: 			return new CellMax();
			case TYPE_MIN: 			return new CellMin();
			case TYPE_KURTOSIS: 	return new CellKurtosis();
			case TYPE_SKEW: 		return new CellSkewness();
			case TYPE_QUANTILE: 	return new CellQuantile(config.quantile, config.quantiles);
			case TYPE_GAP_FRACTION:	return new CellGapFraction(config.gapFractionType);
			default:
				g_argerr("Invalid statistic type: " << config.type);
			}
		}

		void _runner(PointStats *ps) {
			ps->runner();
		}

		void PointStats::runner() {
			size_t idx;
			std::list<std::shared_ptr<LASPoint> > pts;
			while(m_running || m_idxq.size()) {
				if(!m_idxq.size())
					continue;
				m_qmtx.lock();
				if(!m_idxq.size()) {
					m_qmtx.unlock();
					continue;
				}
				idx = m_idxq.front();
				pts = m_cache[idx];
				m_idxq.pop();
				m_qmtx.unlock();
				
				//g_debug(" -- computing stats for " << idx << "; " << pts.size());
				m_mem->set(idx, m_computer->compute(pts));
				m_cmtx.lock();
				m_cache.erase(idx);
				m_cmtx.unlock();
			}
			
		}	
				
		void PointStats::pointstats(const PointStatsConfig &config, const Callbacks *callbacks) {

			checkConfig(config);
			
		 	if(config.threads > 0) {
		 		g_debug(" -- pointstats running with " << config.threads << " threads");
		 		omp_set_dynamic(1);
		 		omp_set_num_threads(config.threads);
		 	} else {
		 		g_argerr("Run with >=1 threads.");
		 	}

			std::vector<std::string> files(config.sourceFiles.begin(), config.sourceFiles.end());
			FinalizedPointStream ps(files, g_abs(config.resolution));
			Bounds bounds = ps.bounds();
			//size_t pointCount = ps.pointCount();
			g_debug(" -- pointstats - work bounds: " << bounds.print() << "; " << ps.cols() << "x" << ps.rows());

			// Prepare the grid
			// TODO: Only works with UTM north.
			Raster<float> grid(config.dstFile, 1, bounds, config.resolution, 
				-config.resolution, -9999, config.hsrid);
			g_debug(" -- pointstats - raster size: " << grid.cols() << ", " << grid.rows());

			m_mem.reset(new MemRaster<float>(ps.cols(), ps.rows(), true));
			m_mem->fill(-9999.0);

			m_computer.reset(getComputer(config));
			CellStatsFilter *filter = nullptr;
			if(config.hasClasses()) 
				filter = new ClassFilter(config.classes);
			//if(config.hasQuantileFilter())
			//	tmp = tmp->chain(new QuantileFilter(config.quantileFilter, config.quantileFilterFrom, config.quantileFilterTo));
			if(filter)
				m_computer->setFilter(filter);

			std::list<std::thread> threads;

			m_running = true;
			size_t finalIdx = 0;

			for(uint32_t i = 0; i < config.threads; ++i) {
				std::thread t(_runner, this);
				threads.push_back(std::move(t));
			}

			g_debug(" -- streaming points");
			LASPoint pt;
			while(ps.next(pt, &finalIdx)) {
				std::shared_ptr<LASPoint> up(new LASPoint(pt));
				m_cmtx.lock();
				m_cache[ps.toIdx(pt)].push_back(up);
				m_cmtx.unlock();
				if(finalIdx) {
					m_qmtx.lock();
					m_idxq.push(finalIdx);
					m_qmtx.unlock();
				}
			}
			
			m_running = false;

			for(std::thread &t : threads)
				t.join();

			if(config.normalize)
				m_mem->normalize();

			g_debug(" -- writing to output");
			grid.writeBlock(*(m_mem.get()));

			if(callbacks)
				callbacks->overallCallback(1.0f);

		}

	} // point

} // geotools

