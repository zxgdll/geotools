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
#include <math.h>
#include <exception>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "geos/index/quadtree/Quadtree.h"

#include <liblas/liblas.hpp>

#include "lasgrid.hpp"
#include "lasutil.hpp"
#include "sortedpointstream.hpp"
#include "simplegeom.hpp"

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;

using namespace geotools::util;
using namespace geotools::raster;
using namespace geotools::las;
using namespace geotools::geom;

namespace geotools {

	namespace las {

		namespace lasgrid_config {
			
			double defaultResolution = 2.0; 
			double defaultRadius = std::sqrt(g_sq(defaultResolution / 2.0) * 2.0); 
			bool defaultSnapToGrid = true;
			unsigned char defaultType = TYPE_MEAN;
			unsigned char defaultAttribute = ATT_HEIGHT;
			unsigned char defaultAngleLimit = 180;
			std::set<unsigned char> defaultClasses = {2}; 
			std::map<std::string, unsigned char> types = {
				{"Minimum", TYPE_MIN}, {"Maximum", TYPE_MAX}, {"Mean", TYPE_MEAN}, {"Density", TYPE_DENSITY},
				{"Sample Variance", TYPE_VARIANCE}, {"Sample Std. Dev.", TYPE_STDDEV}, {"Population Variance", TYPE_PVARIANCE},
				{"Population Std. Dev.", TYPE_PSTDDEV}, {"Count", TYPE_COUNT}, {"Quantile", TYPE_QUANTILE}, {"Median", TYPE_MEDIAN}
			};
			std::map<std::string, unsigned char> attributes = {
				{"Height", ATT_HEIGHT}, {"Intensity", ATT_INTENSITY}
			};

		} // config
		
		namespace lasgrid_util {

			/**
			 * Comparator for sorting doubles.
			 */
			int _fcmp(const void * a, const void * b) {
				const double * aa = (const double *) a;
				const double * bb = (const double *) b;
				return (*aa > *bb) - (*bb > *aa);
			}

			void vector_dealloc(std::vector<double> *item) {
				delete item;
			}

			/**
			 * Returns true if the point is within the radius associated with a cell's centroid.
			 * @param px The x coordinate of the point.
			 * @param py The y coordinate of the point.
			 * @param col The column of the cell of interest.
			 * @param row The row of the cell of interest.
			 * @param radius The radius around the cell's centroid.
			 * @param resolution The resolution of the output raster.
			 * @param bounds The bounds of the raster.
			 */
			bool inRadius(double px, double py, int col, int row, double radius,
					double resolution, Bounds &bounds){
				if(radius == 0.0) return true;
				// If a radius is given, extract the x and y of the current cell's centroid
				// and measure its distance (squared) from the point.
				double x = col * resolution + bounds[0] + resolution * 0.5;
				double y = row * resolution + bounds[1] + resolution * 0.5;
				// If the cell is outside the radius, ignore it.
				double r = sqrt(g_sq(x - px) + g_sq(y - py));
				return r <= radius;
			}

		} // util

		unsigned char LASGridConfig::parseAtt(const std::string &attStr) {
			if("intensity" == attStr) {
				return ATT_INTENSITY;
			} else if("height" == attStr) {
				return ATT_HEIGHT;
			} 
			return 0;
		}

		unsigned char LASGridConfig::parseType(const std::string &typeStr) {
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

		void LASGrid::setCallbacks(Callbacks *callbacks) {
			m_callbacks.reset(callbacks);
		}

		void LASGrid::checkConfig(const LASGridConfig &config) {
			if(config.resolution <= 0.0)
				g_argerr("Resolution must be > 0: " << config.resolution);
			if(config.radius <= 0.0)
				g_argerr("Radius invalid: " << config.radius);
			if(config.lasFiles.size() == 0)
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

			g_debug("Radius: " << config.radius);
			g_debug("Resolution: " << config.resolution);
			g_debug("Files: " << config.lasFiles.size());
			g_debug("Destination: " << config.dstFile);
			g_debug("Attribute: " << config.attribute);
			g_debug("Type: " << config.type);
			g_debug("Classes: " << config.classes.size());
			g_debug("Angle Limit: " << config.angleLimit);
		}

		void LASGrid::computeWorkBounds(const std::list<std::string> &files, const Bounds &bounds,
			std::set<std::string> &selectedFiles, Bounds &workBounds) {
			g_debug("Work bounds initial: " << workBounds.print());
			liblas::ReaderFactory rf;
			for(const std::string &file : files) {
				g_debug("Checking file " << file);
				std::ifstream in(file.c_str(), std::ios::in|std::ios::binary);
				liblas::Reader r = rf.CreateWithStream(in);
				liblas::Header h = r.GetHeader();
				Bounds lasBounds;
				if(!LasUtil::computeLasBounds(h, lasBounds, 2))
					LasUtil::computeLasBounds(r, lasBounds, 2); // If the header bounds are bogus.
				g_debug("File bounds " << file << ": " << lasBounds.print());
				in.close();
				if(bounds.intersects(lasBounds, 2)) {
					selectedFiles.insert(file);
					workBounds.extend(lasBounds);
				}
			}
			g_debug("Work bounds final: " << workBounds.print());
		}


		class Cell {
		public:
			std::unique_ptr<geos::geom::Polygon> circle;
			double x;
			double y;
			double radius;
			std::list<LASPoint> points;

			Cell(double x, double y, double radius) :
				x(x), y(y), radius(radius) {
				circle = std::move(SimpleGeom::createCircle(geos::geom::Coordinate(x, y), radius));
			}

			bool contains(double x, double y) {
				return circle->contains(SimpleGeom::createPoint(x, y, 0.0).get());
			}
		};

		double compute(const Cell *cell) {
			double sum = 0.0;
			int count = cell->points.size();
			for(const LASPoint &pt : cell->points)
				sum += pt.z;
			return count > 0 ? sum / cell->points.size() : 1000.0;
		}

		void LASGrid::lasgrid(const LASGridConfig &config) {

			checkConfig(config);
			
			// Stores the selected files.
			std::set<std::string> files;
			// Stores the bounds of the work area.
			Bounds workBounds;
			workBounds.collapse();

			// Compute the work bounds and indices.
			computeWorkBounds(config.lasFiles, config.bounds, files, workBounds);

			// Snap the bounds if necessary.
			if(config.snap) {
				workBounds.snap(config.resolution);
				g_debug(" -- lasgrid - snapped work bounds: " << workBounds.print());
			}

			// To keep track of completed work.
			geos::geom::Geometry *geom = nullptr;

			// To store cells.
			std::unordered_map<unsigned long, Cell*> cells;

			// Prepare the grid
			// TODO: Only works with UTM north.
			Raster<float> grid(config.dstFile, 1, workBounds, config.resolution, 
				-config.resolution, -9999, config.hsrid);
			g_debug(" -- lasgrid - raster size: " << grid.cols() << ", " << grid.rows());

			geos::index::quadtree::Quadtree tree;

			// Iterate over selected files to produce grid
			float curFile = 0.0f;
			for(const std::string &file : files) {

				g_debug(" -- lasgrid - processing " << file);

				if(m_callbacks)
					m_callbacks->overallCallback((curFile + 0.5f) / files.size());

				SortedPointStream pts(file, 100);
				pts.init();
				
				LASPoint pt;
				unsigned int curPt = 0;
				while(pts.next(pt, &geom)) {

					if(curPt++ % 1000 == 0) {
						if(m_callbacks) 
							m_callbacks->fileCallback((curPt + 1.0f) / pts.pointCount());
					}

					// If the point is outside the scan angle range, skip it.
					if(g_abs(pt.angle) > config.angleLimit || !config.hasClass(pt.cls))
						continue;

					int col = grid.toCol(pt.x); 
					int row = grid.toRow(pt.y); 
					unsigned long idx = ((unsigned long) col << 32) | ((unsigned long) row);

					//g_debug(" -- lasgrid - idx " << idx);
					if(cells.find(idx) == cells.end()) {
						//g_debug(" -- lasgrid - new cell at " << col << "," << row);
						Cell *c = new Cell(
							grid.toX(col) + grid.resolutionX() * 0.5, 
							grid.toY(row) + grid.resolutionY() * 0.5, 
							config.radius
						);
						tree.insert(c->circle->getEnvelopeInternal(), c);
						cells[idx] = c;
					}
					
					// find cells in tree that contain the point
					// add point to cells
					std::vector<void*> res;
					geos::geom::Envelope e(pt.x, pt.y, pt.x, pt.y);
					tree.query(&e, res);

					g_debug(" -- lasgrid - found " << res.size() << " cells for point " << pt.x << "," << pt.y);
					for(const void *c : res)
						((Cell *) c)->points.push_back(pt);

					if(curPt % 1000 == 0 || curPt == pts.pointCount()) {
						g_debug(" -- lasgrid - finalizing cells");
						std::set<unsigned long> rem;
						for(const auto &it : cells) {
							if(geom->contains(it.second->circle.get())) {
								grid.set(grid.toCol(it.second->x), grid.toRow(it.second->y), compute(it.second));
								rem.insert(it.first);
								tree.remove(it.second->circle->getEnvelopeInternal(), it.second);
							}
						}
						g_debug(" -- lasgrid - removing " << rem.size() << " items");
						for(const unsigned long &idx : rem)
							cells.erase(idx);
					}

				}

				if(m_callbacks)
					m_callbacks->overallCallback((curFile + 1.0f) / files.size());
				curFile += 1.0f;
			}

			g_debug(" -- lasgrid - finalizing cells");
			std::set<unsigned long> rem;
			for(const auto &it : cells) {
				if(geom->contains(it.second->circle.get())) {
					grid.set(grid.toCol(it.second->x), grid.toRow(it.second->y), compute(it.second));
					rem.insert(it.first);
					tree.remove(it.second->circle->getEnvelopeInternal(), it.second);
				}
			}
			g_debug(" -- lasgrid - removing " << rem.size() << " items");
			for(const unsigned long &idx : rem)
				cells.erase(idx);

			if(m_callbacks)
				m_callbacks->overallCallback(1.0f);

		}

	} // las

} // geotools

