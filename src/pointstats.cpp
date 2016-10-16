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

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "pointstats.hpp"
#include "lasutil.hpp"
#include "pointstream.hpp"

#include <CGAL/Plane_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2_algorithms.h>

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;

using namespace geotools::util;
using namespace geotools::raster;
using namespace geotools::point;
using namespace geotools::las;
using namespace geotools::point::pointstats_util;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
typedef CGAL::Projection_traits_xy_3<K>  					Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> 					Delaunay;
typedef K::Point_3 											Point_3;
typedef K::Plane_3 											Plane_3;
typedef Delaunay::Finite_faces_iterator						Finite_faces_iterator;
typedef Delaunay::Face										Face;

namespace geotools {

	namespace point {

		namespace pointstats_config {
			
			double defaultResolution = 2.0; 
			bool defaultSnapToGrid = true;
			unsigned char defaultType = TYPE_MEAN;
			unsigned char defaultAttribute = ATT_HEIGHT;
			unsigned char defaultAngleLimit = 180;
			unsigned char defaultGapFraction = GAP_BLB;
			std::set<unsigned char> defaultClasses = {2}; 

			std::map<std::string, unsigned char> types = {
				{"Minimum", TYPE_MIN}, {"Maximum", TYPE_MAX}, {"Mean", TYPE_MEAN}, {"Density", TYPE_DENSITY},
				{"Sample Variance", TYPE_VARIANCE}, {"Sample Std. Dev.", TYPE_STDDEV}, {"Population Variance", TYPE_PVARIANCE},
				{"Population Std. Dev.", TYPE_PSTDDEV}, {"Count", TYPE_COUNT}, {"Quantile", TYPE_QUANTILE}, 
				{"Median", TYPE_MEDIAN}, {"Rugosity", TYPE_RUGOSITY}, {"Kurtosis", TYPE_KURTOSIS}, {"Skewness", TYPE_SKEW},
				{"Gap Fraction", TYPE_GAP_FRACTION}
			};
			std::map<std::string, unsigned char> attributes = {
				{"Height", ATT_HEIGHT}, {"Intensity", ATT_INTENSITY}
			};
			std::map<std::string, unsigned char> gapFractionTypes = {
				{"IR", GAP_IR}, {"BLa", GAP_BLA}, {"BLb", GAP_BLB}, {"RR", GAP_RR}, {"FR", GAP_FR}
			};

		} // config
		
		namespace pointstats_util {

			CellStats::~CellStats() {}

			class CellDensity : public CellStats {
			private:
				double m_cellArea;
			public:
				CellDensity(double cellArea) :
					m_cellArea(cellArea) {
				}

				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					if(values.size() == 0)
						return -9999.0;
					return values.size() / m_cellArea;
				}
			};

			class CellMean : public CellStats {
			public:
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					if(values.size() == 0)
						return -9999.0;
					double sum = 0.0;
					for(const std::unique_ptr<LASPoint> &v : values)
						sum += v->z;
					return sum / values.size();
				}
			};

			class CellCount : public CellStats {
			public:
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					return values.size();
				}
			};

			class CellMedian : public CellStats {
			public:
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					if(values.size() <= 1)
						return -9999.0;
					int i = 0;
					std::vector<double> v(values.size());
					for(const std::unique_ptr<LASPoint> &pt : values)
						v[i++] = pt->z;
					std::sort(v.begin(), v.end());
					unsigned int size = v.size();
					if(size % 2 == 0) {
						return (v[(int) size / 2] + v[(int) size / 2 - 1]) / 2.0;
					} else {
						return v[(int) size / 2];
					}
				}
			};

			class CellMin : public CellStats {
			public:
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					if(values.size() <= 1)
						return -9999.0;
					double min = G_DBL_MAX_POS;
					for(const std::unique_ptr<LASPoint> &v : values) {
						if(v->z < min)
							min = v->z;
					}
					return min;
				}
			};

			class CellMax : public CellStats {
			public:
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					if(values.size() <= 1)
						return -9999.0;
					double max = G_DBL_MAX_NEG;
					for(const std::unique_ptr<LASPoint> &v : values) {
						if(v->z > max)
							max = v->z;
					}
					return max;
				}
			};

			class CellSampleVariance : public CellStats {
			private:
				CellMean m_mean;
			public:
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					if(values.size() <= 1)
						return -9999.0;			
					double mean = m_mean.compute(values);
					double sum = 0;
					for(const std::unique_ptr<LASPoint> &v : values)
						sum += g_sq(g_abs(v->z - mean));
					return sum / (values.size() - 1);
				}
			};

			class CellPopulationVariance : public CellStats {
			private:
				CellMean m_mean;
			public:
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					if(values.size() <= 1)
						return -9999.0;			
					double mean = m_mean.compute(values);
					double sum = 0;
					for(const std::unique_ptr<LASPoint> &v : values)
						sum += g_sq(g_abs(v->z - mean));
					return sum / values.size();
				}
			};

			class CellSampleStdDev : public CellStats {
			private:
				CellSampleVariance m_variance;
			public:
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					if(values.size() == 0)
						return -9999.0;
					return std::sqrt(m_variance.compute(values));
				}
			};

			class CellPopulationStdDev : public CellStats {
			private:
				CellPopulationVariance m_variance;
			public:
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					if(values.size() == 0)
						return -9999.0;
					return std::sqrt(m_variance.compute(values));
				}
			};

			class CellSkewness : public CellStats {
			private:
				CellMean m_mean;
				CellSampleStdDev m_stdDev;
			public:
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					if(values.size() == 0)
						return -9999.0;
					// Fisher-Pearson
					double mean = m_mean.compute(values);
					double sum = 0.0;
					unsigned int count = values.size();
					for(const std::unique_ptr<LASPoint> &v : values)
						sum += std::pow(v->z - mean, 3.0) / count;
					return sum / std::pow(m_stdDev.compute(values), 3.0);
				}
			};

			class CellKurtosis : public CellStats {
			private:
				CellMean m_mean;
				CellSampleStdDev m_stdDev;
			public:
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					if(values.size() == 0)
						return -9999.0;
					double mean = m_mean.compute(values);
					double sum = 0.0;
					unsigned int count = values.size();
					for(const std::unique_ptr<LASPoint> &v : values)
						sum += std::pow(v->z - mean, 4.0) / count;
					return sum / std::pow(m_stdDev.compute(values), 4.0) - 3.0;
				}
			};

			class CellQuantile : public CellStats {
			private:
				unsigned int m_quantile;
				unsigned int m_quantiles;
			public:
				CellQuantile(unsigned char quantile, unsigned char quantiles) :
					m_quantile(quantile), m_quantiles(quantiles) {
				}
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					return 0;
				}
			};



			class CellRugosity : public CellStats {
			private:
				static const int a = 2, b = 1, c = 0;

				double computeArea(const Face &face) {
					Point_3 p1 = face.vertex(0)->point();
					Point_3 p2 = face.vertex(1)->point();
					Point_3 p3 = face.vertex(2)->point();
					double sides[] {
						std::sqrt(std::pow(p1.x() - p2.x(), 2.0) + std::pow(p1.y() - p2.y(), 2.0) + std::pow(p1.z() - p2.z(), 2.0)),
						std::sqrt(std::pow(p2.x() - p3.x(), 2.0) + std::pow(p2.y() - p3.y(), 2.0) + std::pow(p2.z() - p3.z(), 2.0)),
						std::sqrt(std::pow(p3.x() - p1.x(), 2.0) + std::pow(p3.y() - p1.y(), 2.0) + std::pow(p3.z() - p1.z(), 2.0)),
					};
					double s = (sides[a] + sides[b] + sides[c]) / 2.0;
					return std::sqrt(s * (s - sides[a]) * (s - sides[b]) * (s - sides[c]));
				}

				double toPlane(const Point_3 &p, const Plane_3 &plane, const Point_3 &centroid) {
					return (p.x() * plane.a() + p.y() * plane.b() + plane.d()) / -plane.c();
				}

			public:

				/**
				 * Using Du Preez, 2014 - Arc-Chord Ratio (ACR) Index.
				 */
				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {

					if(values.size() == 0)
						return -9999.0;
					
					std::list<Point_3> pts;
					for(const std::unique_ptr<LASPoint> &v : values)
						pts.push_back(Point_3(v->x, v->y, v->z));
					
					// Delaunay 3D surface area.
					double tarea = 0.0;
					Delaunay dt(pts.begin(), pts.end());
					for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it)
						tarea += computeArea(*it);

					// Convex hull and POBF
					std::list<Point_3> hull;
					Plane_3 plane;
					Point_3 centroid;
					CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(hull), Gt());
					CGAL::linear_least_squares_fitting_3(hull.begin(), hull.end(), plane, centroid, CGAL::Dimension_tag<0>());

					// POBF surface area.
					double parea = 0.0;
					std::list<Point_3> poly; // TODO: Is there a faster way?
					for(const Point_3 &pt : hull)
						poly.push_back(Point_3(pt.x(), pt.y(), toPlane(pt, plane, centroid)));
					Delaunay dt2(poly.begin(), poly.end());
					for(Finite_faces_iterator it = dt2.finite_faces_begin(); it != dt2.finite_faces_end(); ++it)
						parea += computeArea(*it);

					return tarea / parea;
				}
			};

			/**
			 * Adapted from:
			 * Hopkinson, C., & Chasmer, L. (2009). Testing LiDAR models of fractional cover across 
			 *	multiple forest ecozones. Remote Sensing of Environment, 113(1), 275–288. 
			 *	http://doi.org/10.1016/j.rse.2008.09.012
			 */
			class CellGapFraction : public CellStats {
			private:
				unsigned char m_type;

				double fcLidarBLa(const std::list<std::unique_ptr<LASPoint> > &values) {
					double gnd = 0.0, all = 0.0;
					for(const std::unique_ptr<LASPoint> &pt : values) {
						if(pt->ground())
							gnd += pt->intensity;
						all += pt->intensity; // TODO: This should perhaps be filtered by class to remove bogus points.
					}
					g_debug(" -- fcLidarBLa " << gnd << ", " << all << ", " << values.size());
					return all != 0.0 ? 1.0 - std::sqrt(gnd / all) : -9999.0;
				}

				double fcLidarBLb(const std::list<std::unique_ptr<LASPoint> > &values) {
					double gndSingle = 0.0, gndLast = 0.0, first = 0.0, single = 0.0, intermediate = 0.0, last = 0.0, total = 0.0;
					for(const std::unique_ptr<LASPoint> &pt : values) {
						if(pt->ground()) {
							if(pt->single())
								gndSingle += pt->intensity;
							if(pt->last())
								gndLast += pt->intensity;
						}
						if(pt->first())
							first += pt->intensity;
						if(pt->single())
							single += pt->intensity;
						if(pt->intermediate())
							intermediate += pt->intensity;
						if(pt->last())
							last += pt->intensity;
						total += pt->intensity; // TODO: This should perhaps be filtered by class to remove bogus points.
					}
					if(total == 0.0) return -9999.0;
					double denom = (first + single) / total + std::sqrt((intermediate + last) / total);
					if(denom == 0.0) return -9999.;
					return (gndSingle / total + std::sqrt(gndLast / total)) / denom;
				}

				double fcLidarIR(const std::list<std::unique_ptr<LASPoint> > &values) {
					double canopy = 0.0, total = 0.0;
					for(const std::unique_ptr<LASPoint> &pt : values) {
						if(!pt->ground())
							canopy += pt->intensity;
						total += pt->intensity;
					}
					return total != 0.0 ? canopy / total : -9999.0;
				}

				
				double fcLidarRR(const std::list<std::unique_ptr<LASPoint> > &values) {
					unsigned int canopy = 0, total = 0;
					for(const std::unique_ptr<LASPoint> &pt : values) {
						if(!pt->ground())
							++canopy;
						++total;
					}
					return total != 0.0 ? (double) canopy / total : -9999.0;
				}
				
				double fcLidarFR(const std::list<std::unique_ptr<LASPoint> > &values) {
					unsigned int canopy = 0, total = 0;
					for(const std::unique_ptr<LASPoint> &pt : values) {
						if(pt->first()) {
							if(!pt->ground())
								++canopy;
							++total;
						}
					}
					return total != 0.0 ? (double) canopy / total : -9999.0;
				}
			public:
				const static unsigned char IR = GAP_IR;
				const static unsigned char BLA = GAP_BLA;
				const static unsigned char BLB = GAP_BLB;
				const static unsigned char RR = GAP_RR;
				const static unsigned char FR = GAP_FR;

				CellGapFraction(unsigned char type) :
					m_type(type) {
				}

				double compute(const std::list<std::unique_ptr<LASPoint> > &values) {
					switch(m_type) {
					case BLA:	return fcLidarBLa(values);
					case BLB:	return fcLidarBLb(values);
					case IR:	return fcLidarIR(values);
					case RR:	return fcLidarRR(values);
					case FR:	return fcLidarFR(values);
					default:
						g_argerr("Unknown Gap Fraction method: " << m_type);
					}
				}
			};

		} // util

		unsigned char PointStatsConfig::parseAtt(const std::string &attStr) {
			if("intensity" == attStr) {
				return ATT_INTENSITY;
			} else if("height" == attStr) {
				return ATT_HEIGHT;
			} 
			return 0;
		}

		unsigned char PointStatsConfig::parseType(const std::string &typeStr) {
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

		void PointStats::computeWorkBounds(const std::list<std::string> &files, const Bounds &bounds,
			std::set<std::string> &selectedFiles, Bounds &workBounds, unsigned long *pointCount) {
			g_debug(" -- computeWorkBounds - work bounds initial: " << workBounds.print());
			liblas::ReaderFactory rf;
			unsigned long count = 0;
			for(const std::string &file : files) {
				g_debug(" -- computeWorkBounds - checking file " << file);
				geotools::las::PointStream ps(file);
				count += ps.pointCount();
				if(bounds.intersects(ps.fileBounds(), 2)) {
					selectedFiles.insert(file);
					workBounds.extend(ps.fileBounds());
				}
			}
			*pointCount = count;
			g_debug(" -- computeWorkBounds - work bounds final: " << workBounds.print() << "; point count " << *pointCount);
		}

		CellStats* PointStats::getComputer(const PointStatsConfig &config) {
			switch(config.type) {
			case TYPE_MEAN: 		return new CellMean();
			case TYPE_MEDIAN: 		return new CellMedian();
			case TYPE_COUNT: 		return new CellCount();
			case TYPE_STDDEV: 		return new CellSampleStdDev();
			case TYPE_VARIANCE: 	return new CellSampleVariance();
			case TYPE_PSTDDEV: 		return new CellPopulationStdDev();
			case TYPE_PVARIANCE:	return new CellPopulationVariance();
			case TYPE_DENSITY: 		return new CellDensity(g_sq(config.resolution));
			case TYPE_RUGOSITY: 	return new CellRugosity();
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

		void PointStats::pointstats(const PointStatsConfig &config, const Callbacks *callbacks) {

			checkConfig(config);
			
		 	if(config.threads > 0) {
		 		g_debug(" -- pointstats running with " << config.threads << " threads");
		 		omp_set_dynamic(1);
		 		omp_set_num_threads(config.threads);
		 	} else {
		 		g_argerr("Run with >=1 threads.");
		 	}

			// Compute the work bounds, and store the list
			// of relevant files. Snap the bounds if necessary.
			std::set<std::string> files;
			unsigned long pointCount;
			Bounds workBounds;
			workBounds.collapse();
			computeWorkBounds(config.sourceFiles, config.bounds, files, workBounds, &pointCount);
			if(config.snap) {
				workBounds.snap(config.resolution);
				g_debug(" -- pointstats - snapped work bounds: " << workBounds.print());
			}

			// Prepare the grid
			// TODO: Only works with UTM north.
			Raster<float> grid(config.dstFile, 1, workBounds, config.resolution, 
				-config.resolution, -9999, config.hsrid);
			grid.fill(-9999.0);
			g_debug(" -- pointstats - raster size: " << grid.cols() << ", " << grid.rows());

			std::unique_ptr<CellStats> computer(getComputer(config));
			std::vector<std::string> filesv(files.begin(), files.end());
			unsigned int fileIdx = 0;

			std::map<unsigned long, std::list<std::unique_ptr<LASPoint> > > edges;

			ClassFilter filter(config.classes);

			#pragma omp parallel for
			for(unsigned int f = 0; f < filesv.size(); ++f) {

				g_debug(" -- pointstats - file " << filesv[f]);

				if(callbacks)
					callbacks->overallCallback((fileIdx + 0.5f) / filesv.size());

				const std::string &file = filesv[f];
				std::map<unsigned long, std::list<std::unique_ptr<LASPoint> > > values;

				LASPoint pt;
				PointStream ps(file);

				while(ps.next(pt, &filter)) {
					unsigned long idx = grid.toRow(pt.y) * grid.cols() + grid.toCol(pt.x);
					std::unique_ptr<LASPoint> pv(new LASPoint(pt));
					values[idx].push_back(std::move(pv));
				}

				g_debug(" -- pointstats - computing results");
				std::list<unsigned long> edgeIds;
				for(const auto &it : values) {
					int col = it.first % grid.cols();
					int row = it.first / grid.cols();
					if(ps.contains(grid.toX(col), grid.toY(row)) 
							&& ps.contains(grid.toX(col + 1), grid.toY(row + 1))) {
						grid.set(it.first, computer->compute(it.second));
					} else {
						edgeIds.push_back(it.first);
					}
				}

				g_debug(" -- pointstats - " << edgeIds.size() << " edge cells deferred.");
				// Write to the edge pixel map to solve later.
				#pragma omp critical
				{
					for(const unsigned long &id : edgeIds) {
						for(std::unique_ptr<LASPoint> &pt : values[id])
							edges[id].push_back(std::move(pt));
					}
	
					++fileIdx; // Do it here to save another critical/atomic
				}

				if(callbacks)
					callbacks->overallCallback((float) fileIdx / filesv.size()); // N.B. fileIdx already incremented

			}

			unsigned int i = 0;
			std::vector<unsigned long> ids(edges.size());
			for(const auto &it : edges)
				ids[i++] = it.first;
			std::sort(ids.begin(), ids.end());

			#pragma omp parallel for
			for(unsigned int i = 0; i < ids.size(); ++i) {
				grid.set(ids[i], computer->compute(edges[ids[i]]));
			}

			if(callbacks)
				callbacks->overallCallback(1.0f);

		}

	} // point

} // geotools

