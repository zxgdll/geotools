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

		class Tree {
		public:
			double tx;
			double ty;
			double radius;
			double minRadius;
			std::list<Point> points;

			Tree *tl;
			Tree *tr;
			Tree *bl;
			Tree *br;

			Tree(double x, double y, double radius, double minRadius) :
				tx(x), ty(y), radius(radius), minRadius(minRadius),
				tl(nullptr), tr(nullptr), bl(nullptr), br(nullptr) {

			}

			void getContained(Bounds &cbounds, std::list<Tree*> &trees) {
				if(radius <= minRadius) {
					if(cbounds.contains(tx - radius, ty - radius)
						&& cbounds.contains(tx + radius, ty + radius)
						&& points.size() > 0) {
						trees.push_back(this);
					}
				} else {
					if(tl != nullptr) tl->getContained(cbounds, trees);
					if(bl != nullptr) bl->getContained(cbounds, trees);
					if(tr != nullptr) tr->getContained(cbounds, trees);
					if(br != nullptr) br->getContained(cbounds, trees);
				}
			}

			bool contains(double x, double y) {
				return std::sqrt(g_sq(x - tx) + g_sq(y - ty)) <= radius;
			}

			void getContains(double x, double y, std::list<Tree*> &trees) {
				if(contains(x, y)) {
					if(radius <= minRadius) {
						trees.push_back(this);
					} else {
						if(tl != nullptr) tl->getContains(x, y, trees);
						if(tr != nullptr) tr->getContains(x, y, trees);
						if(bl != nullptr) bl->getContains(x, y, trees);
						if(br != nullptr) br->getContains(x, y, trees);
					}
				}
			}

			void add(double x, double y, double z) {
				if(contains(x, y)) {
					if(radius <= minRadius) {
						g_debug(" -- add point: " << x << "," << y << "," << z);
						points.push_back(Point(x, y, z));
					} else {
						if(tl == nullptr) {
							tl = new Tree(std::cos(G_PI * 0.375) * radius / 2.0, 
								std::sin(G_PI * 0.375) * radius / 2.0, radius / 2.0, minRadius);
						}
						tl->add(x, y, z);
						if(bl == nullptr) {
							bl = new Tree(std::cos(G_PI * 0.625) * radius / 2.0, 
								std::sin(G_PI * 0.625) * radius / 2.0, radius / 2.0, minRadius);
						}
						bl->add(x, y, z);
						if(tr == nullptr) {
							tr = new Tree(std::cos(G_PI * 0.125) * radius / 2.0, 
								std::sin(G_PI * 0.125) * radius / 2.0, radius / 2.0, minRadius);
						}
						tr->add(x, y, z);
						if(br == nullptr) {
							br = new Tree(std::cos(G_PI * 0.875) * radius / 2.0, 
								std::sin(G_PI * 0.875) * radius / 2.0, radius / 2.0, minRadius);
						}
						br->add(x, y, z);
					}
				}
			}

			~Tree() {
				if(tl) delete tl;
				if(tr) delete tr;
				if(bl) delete bl;
				if(br) delete br;
			}

		};

		class Cell {
		public:
			std::unique_ptr<geos::geom::Polygon> circle;
			double x;
			double y;
			double radius;
			std::list<Point> points;

			Cell(double x, double y, double radius) :
				x(x), y(y), radius(radius) {
				circle = std::move(SimpleGeom::createCircle(geos::geom::Coordinate(x, y), radius));
			}

			bool contains(double x, double y) {
				return circle->contains(SimpleGeom::createPoint(x, y, 0.0).get());
			}
		};

		double compute(const Cell *cell) {
			return 0.0;
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
				g_debug("Snapped work bounds: " << workBounds.print());
			}

			// To keep track of completed work.
			std::unique_ptr<geos::geom::Geometry> geom;

			// To store cells.
			std::unordered_map<unsigned long, std::unique_ptr<Cell> > cells;

			// Prepare the grid
			// TODO: Only works with UTM north.
			Raster<float> grid(config.dstFile, 1, workBounds, config.resolution, 
				-config.resolution, -9999, config.hsrid);
			g_debug("Raster size: " << grid.cols() << ", " << grid.rows());

			// Iterate over selected files to produce grid
			float idx = 0.0f;
			for(const std::string &file : files) {

				if(m_callbacks)
					m_callbacks->overallCallback((idx + 0.5f) / files.size());

				SortedPointStream pts(file, 100, SortedPointStream::xyz);
				Bounds lasBounds = pts.fileBounds();
				
				if(geom.get() == nullptr)
					geom = std::move(SimpleGeom::createPolygon(lasBounds));

				unsigned int curPt = 0;
				unsigned int numPts = pts.pointCount();

				Point pt;
				geos::geom::Geometry *g = geom.get();

				while(pts.next(pt, &g)) {

					if(curPt % 1000 == 0) {
						if(m_callbacks) 
							m_callbacks->fileCallback((curPt + 1.0f) / numPts);
					}
					++curPt;

					// If the point is outside the scan angle range, skip it.
					//if(g_abs(pt.GetScanAngleRank()) > config.angleLimit)
					//	continue;

					// If this point is not in the class list, skip it.
					//unsigned char cls = pt.GetClassification().GetClass();
					//if(config.classes.find(cls) != config.classes.end())
					//	continue;

					int col = (int) (((pt.x - grid.minx()) / grid.width()) * grid.resolutionX());
					int row = (int) (((pt.y - grid.miny()) / grid.height()) * grid.resolutionY());
					unsigned long idx = ((unsigned long) col << 32) | ((unsigned long) row);

					if(cells.find(idx) == cells.end()) {
						std::unique_ptr<Cell> c(new Cell(grid.toX(col) + grid.resolutionX() * 0.5, grid.toY(row) + grid.resolutionY() * 0.5, config.radius));
						cells[idx] = std::move(c);
					}
					
					if(cells[idx]->contains(pt.x, pt.y))
						cells[idx]->points.push_back(pt);

				}

				std::set<unsigned long> rem;
				for(const auto &it : cells) {
					if(g->contains(it.second->circle.get())) {
						grid.set(grid.toCol(it.second->x), grid.toRow(it.second->y), compute(it.second.get()));
						rem.insert(it.first);
					}
				}
				for(const unsigned long &idx : rem)
					cells.erase(idx);

				if(m_callbacks)
					m_callbacks->overallCallback((idx + 1.0f) / files.size());

				idx += 1.0f;
			}

			if(m_callbacks)
				m_callbacks->overallCallback(1.0f);

		}

	} // las

} // geotools

