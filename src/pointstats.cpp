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

		std::vector<uint8_t> PointStatsConfig::parseTypes(const std::vector<std::string> &typeStrs) {
			std::vector<uint8_t> types;
			for(const std::string &typeStr : typeStrs) {
				if("min" == typeStr) {
					types.push_back(TYPE_MIN);
				} else if("max" == typeStr) {
					types.push_back(TYPE_MAX);
				} else if("mean" == typeStr) {
					types.push_back(TYPE_MEAN);
				} else if("density" == typeStr) {
					types.push_back(TYPE_DENSITY);
				} else if("variance" == typeStr) {
					types.push_back(TYPE_VARIANCE);
				} else if("stddev" == typeStr) {
					types.push_back(TYPE_STDDEV);
				} else if("pvariance" == typeStr) {
					types.push_back(TYPE_PVARIANCE);
				} else if("pstddev" == typeStr) {
					types.push_back(TYPE_PSTDDEV);
				} else if("count" == typeStr) {
					types.push_back(TYPE_COUNT);
				} else if("median" == typeStr) {
					types.push_back(TYPE_MEDIAN);
				} else if("skew" == typeStr) {
					types.push_back(TYPE_SKEW);
				} else if("rugosity" == typeStr) {
					types.push_back(TYPE_RUGOSITY);
				} else if("kurtosis" == typeStr) {
					types.push_back(TYPE_KURTOSIS);
				} else if("gap" == typeStr) {
					types.push_back(TYPE_GAP_FRACTION);
				}
			}
			return types;
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
			if(!config.sourceFiles.size())
				g_argerr("At least one input file is required.");
			if(!config.dstFiles.size()) 
				g_argerr("At least one output file is required.");
			if(config.attribute == 0)
				g_argerr("An attribute is required.");
			if(!config.types.size())
				g_argerr("At least one valid type is required.");
			if(config.classes.size() == 0)
				g_warn("No classes given. Matching all classes.");
			if(config.angleLimit <= 0)
				g_argerr("Angle limit must be greater than zero.");
			if(config.dstFiles.size() != config.types.size())
				g_argerr("There should be one output file for each type.");

			g_debug("Resolution: " << config.resolution);
			g_debug("Files: " << config.sourceFiles.size());
			g_debug("Destinations: " << config.dstFiles.size());
			g_debug("Attribute: " << config.attribute);
			g_debug("Types: " << config.types.size());
			g_debug("Classes: " << config.classes.size());
			g_debug("Angle Limit: " << config.angleLimit);
		}

		CellStats* PointStats::getComputer(const uint8_t &type, const PointStatsConfig &config) {
			using namespace geotools::point::stats;
			switch(type) {
			case TYPE_MEAN: 		return new CellMean();
			case TYPE_MEDIAN: 		return new CellMedian();
			case TYPE_COUNT: 		return new CellCount();
			case TYPE_STDDEV: 		return new CellSampleStdDev();
			case TYPE_VARIANCE: 	return new CellSampleVariance();
			case TYPE_PSTDDEV: 		return new CellPopulationStdDev();
			case TYPE_PVARIANCE:	return new CellPopulationVariance();
			case TYPE_DENSITY: 		return new CellDensity(g_sq(config.resolution));
			case TYPE_RUGOSITY: 	return new CellRugosity(g_sq(config.resolution), 20.0);
			case TYPE_MAX: 			return new CellMax();
			case TYPE_MIN: 			return new CellMin();
			case TYPE_KURTOSIS: 	return new CellKurtosis();
			case TYPE_SKEW: 		return new CellSkewness();
			case TYPE_QUANTILE: 	return new CellQuantile(config.quantile, config.quantiles);
			case TYPE_GAP_FRACTION:	return new CellGapFraction(config.gapFractionType);
			default:
				g_argerr("Invalid statistic type: " << type);
			}
		}

		void _runner(PointStats *ps) {
			ps->runner();
		}

		void PointStats::runner() {
			size_t idx;
			std::list<LASPoint> pts;
			while(m_running) {
				if(!m_bq.pop(&idx))
					continue;
				{
					std::unique_lock<std::mutex> lk(m_cmtx);
					pts.assign(m_cache[idx].begin(), m_cache[idx].end());
					m_cache.erase(idx);
				}
				if(!pts.empty()) {
					for(size_t i = 0; i < m_computers.size(); ++i) {
						std::unique_lock<std::mutex> flk(*(m_mtx[i].get()));
						m_mem[i]->set(idx, m_computers[i]->compute(pts));
					}
					pts.clear();
				}
			}			
		}	
			
		void PointStats::pointstats(const PointStatsConfig &config, const Callbacks *callbacks) {

			checkConfig(config);
			
		 	if(config.threads > 0) {
		 		g_debug(" -- pointstats running with " << g_max(1, config.threads-1) << " threads");
		 	} else {
		 		g_argerr("Run with >=1 threads.");
		 	}

			// Initialize the point stream; it expects a vector. 
			std::vector<std::string> files(config.sourceFiles.begin(), config.sourceFiles.end());
			FinalizedPointStream ps(files, g_abs(config.resolution));
			Bounds bounds = ps.bounds();
			g_debug(" -- pointstats - work bounds: " << bounds.print() << "; " << ps.cols() << "x" << ps.rows());

			// Initialize the filter group.
			CellStatsFilter *filter = nullptr;
			if(config.hasClasses()) 
				filter = new ClassFilter(config.classes);
			//if(config.hasQuantileFilter())
			//	tmp = tmp->chain(new QuantileFilter(config.quantileFilter, config.quantileFilterFrom, config.quantileFilterTo));

			// Create a computer, grid and mutex for each statistic.
			for(size_t i = 0; i < config.types.size(); ++i) {
				// Create computers for each stat
				g_debug(" -- configuring computer " << (int) config.types[i]);
				std::unique_ptr<CellStats> cs(getComputer(config.types[i], config));
				if(filter)
					cs->setFilter(filter);
				m_computers.push_back(std::move(cs));
				
				// Create raster grid for each stat.
				g_debug(" -- configuring grid");
				std::unique_ptr<MemRaster<float> > mr(new MemRaster<float>(ps.cols(), ps.rows(), true));
				mr->fill(-9999.0);
				mr->nodata(-9999.0);
				m_mem.push_back(std::move(mr));

				// Create mutex for grid.
				g_debug(" -- creating mutex");
				std::unique_ptr<std::mutex> m(new std::mutex());
				m_mtx.push_back(std::move(m));
			}

			// Initialize the thread group for runner.
			std::list<std::thread> threads;
			m_running = true;
			size_t finalIdx = 0;

			// Start the runner threads.
			for(uint32_t i = 0; i < g_max(1, config.threads - 1); ++i) {
				std::thread t(_runner, this);
				threads.push_back(std::move(t));
			}

			// Begin streaming the points into the cache for processing.
			g_debug(" -- streaming points");
			LASPoint pt;
			while(ps.next(pt, &finalIdx)) {
				{
					std::unique_lock<std::mutex> lk(m_cmtx);
					m_cache[ps.toIdx(pt)].push_back(pt);
				}
				if(finalIdx)
					m_bq.push(finalIdx);
			}

			m_bq.flush();
			m_running = false;
			m_bq.finish();
			g_debug("done");

			// Shut down and join the runners.
			for(std::thread &t : threads)
				t.join();

			// Write the grids to files.
			for(size_t i = 0; i < m_mem.size(); ++i) {
				// Normalize the grids if desired.
				if(config.normalize) {
					g_debug(" -- normalizing");
					m_mem[i]->normalize();
				}
				// Prepare the grid
				// TODO: Only works with UTM north.
				g_debug(" -- preparing " << config.dstFiles[i]);
				Raster<float> grid(config.dstFiles[i], 1, bounds, config.resolution, 
					-config.resolution, -9999, config.hsrid);

				// Write grid to file.
				g_debug(" -- writing " << config.types[i]);
				grid.writeBlock(*(m_mem[i].get()));
			}

			if(callbacks)
				callbacks->overallCallback(1.0f);

		}

	} // point

} // geotools

