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

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

/*
#include <CGAL/Plane_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2_algorithms.h>
*/

#include "pointstats.hpp"
#include "lasutil.hpp"
#include "sortedpointstream.hpp"
#include "cellstats.hpp"

using namespace geotools::util;
using namespace geotools::raster;
using namespace geotools::point;
using namespace geotools::las;
using namespace geotools::point::stats;

/*
typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
typedef CGAL::Projection_traits_xy_3<K>  					Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> 					Delaunay;
typedef K::Point_3 											Point_3;
typedef K::Plane_3 											Plane_3;
typedef Delaunay::Finite_faces_iterator						Finite_faces_iterator;
typedef Delaunay::Face										Face;
*/

namespace geotools {

	namespace point {

		namespace pointstats_config {
			
			double defaultResolution = 10.0; 
			bool defaultSnapToGrid = true;
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

		void _normalize(MemRaster<float> &mem) {
			double sum = 0.0;
			for(uint64_t i = 0; i < mem.size(); ++i) {
				double v = mem[i];
				if(v != -9999.0 && !std::isnan(v))
					sum += mem[i];
			}
			double mean = sum / mem.size();
			sum = 0.0;
			for(uint64_t i = 0; i < mem.size(); ++i) {
				double v = mem[i];
				if(v != -9999.0 && !std::isnan(v))
					sum += std::pow(mem[i] - mean, 2.0);
			}
			double stdDev = std::sqrt(sum);
			for(uint64_t i = 0; i < mem.size(); ++i) {
				double v = mem[i];
				if(v != -9999.0 && !std::isnan(v))
					mem.set(i, (mem[i] - mean) / stdDev);
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

			SortedPointStream ps(config.sourceFiles, "cache.tmp", -config.resolution * 2.0, config.rebuild);
			ps.init();
			Bounds workBounds = ps.bounds();
			g_debug(" -- pointstats - work bounds: " << workBounds.print());

			// Prepare the grid
			// TODO: Only works with UTM north.
			Raster<float> grid(config.dstFile, 1, workBounds, config.resolution, 
				-config.resolution, -9999, config.hsrid);
			g_debug(" -- pointstats - raster size: " << grid.cols() << ", " << grid.rows());

			MemRaster<float> mem(grid.cols(), grid.rows(), true);
			mem.fill(-9999.0);

			#pragma omp parallel
			{

				std::unique_ptr<CellStats> computer(getComputer(config));
				CellStatsFilter *filter = nullptr, *tmp;
				if(config.hasClasses()) 
					tmp = filter = new ClassFilter(config.classes);
				//if(config.hasQuantileFilter())
				//	tmp = tmp->chain(new QuantileFilter(config.quantileFilter, config.quantileFilterFrom, config.quantileFilterTo));
				if(filter)
					computer->setFilter(filter);

				#pragma omp for
				for(uint32_t i = 0; i < ps.rowCount(); ++i) {

					//g_debug(" -- row " << i << "; thread " << omp_get_thread_num());
					std::list<std::shared_ptr<LASPoint> > row;
					#pragma omp critical(__row)
					ps.next(row);

					if(row.size()) {

						std::unordered_map<uint64_t, std::list<std::shared_ptr<LASPoint> > > values;

						for(const std::shared_ptr<LASPoint> &pt : row) {
							uint64_t idx = grid.toRow(pt->y) * grid.cols() + grid.toCol(pt->x);
							//g_debug(" -- " << grid.toCol(pt->x) << ", " << grid.toRow(pt->y));
							values[idx].push_back(pt);
						}

						for(const auto &it : values) {
							//g_debug(" -- index " << it.first);
							float val = computer->compute(it.second);
							uint64_t idx = it.first;
							uint32_t col = idx % grid.cols();
							uint32_t row = idx / grid.cols();
							//g_debug(" -- " << col << ", " << row);
							mem.set(it.first, val);
						}

					}
				}
			}

			//_normalize(mem);

			g_debug(" -- writing to output");
			grid.writeBlock(mem);

			if(callbacks)
				callbacks->overallCallback(1.0f);

		}

	} // point

} // geotools

