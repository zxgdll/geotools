#ifndef __SORTEDPOINTSTREAM_HPP__
#define __SORTEDPOINTSTREAM_HPP__

#include <map>
#include <string>
#include <fstream>
#include <memory>
#include <unordered_map>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Boolean_set_operations_2.h>

#include "liblas/liblas.hpp"

#include "geos/geom/Geometry.h"

#include "util.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 											Point_2;
typedef CGAL::Polygon_with_holes_2<K> 						Polygon_h_2;
typedef CGAL::Polygon_2<K>	 								Polygon_2;

using namespace geotools::util;

namespace geotools {
	namespace las {

		class BoundsTracker {
		private:			
			Polygon_h_2 *m_bounds;

		public:

			BoundsTracker() : 
				m_bounds(nullptr) {}

			void add(double x1, double y1, double x2, double y2) {
				if(x1 == x2 || y1 == y2)
					g_argerr("X and Y coords must differ.");
				if(x1 > x2) {
					double tmp = x1;
					x1 = x2;
					x2 = tmp;
				}
				if(y1 > y2) {
					double tmp = y1;
					y1 = y2;
					y2 = tmp;
				}
				std::list<Point_2> pts;
				pts.push_back(Point_2(x1, y1));
				pts.push_back(Point_2(x1, y2));
				pts.push_back(Point_2(x2, y2));
				pts.push_back(Point_2(x2, y1));
				pts.push_back(Point_2(x1, y1));
				if(!m_bounds) {
					m_bounds = new Polygon_h_2(Polygon_2(pts.begin(), pts.end()));
				} else {
					Polygon_2 input(pts.begin(), pts.end());
					Polygon_h_2 *output = new Polygon_h_2();
					CGAL::intersection(*m_bounds, input, output);
					delete m_bounds;
					m_bounds = output;
				}
			}

			bool contains(double x, double y) const {
				return CGAL::ON_BOUNDED_SIDE == m_bounds->outer_boundary().bounded_side(Point_2(x, y));
			}

			bool contains(double x1, double y1, double x2, double y2) const {
				return contains(x1, y1) && contains(x2, y2);
			}

			~BoundsTracker() {
				delete m_bounds;
			}
		};

		class LASPoint {
		public:
			static double scaleX, scaleY, scaleZ;
			double x, y, z;
			uint16_t intensity;
			uint16_t ret, numRets;
			uint16_t scanDir;
			uint8_t cls;
			int8_t angle;

			LASPoint() {}
			LASPoint(const liblas::Point &pt) {
				update(pt);
			}

			static void setScale(double x, double y, double z) {
				scaleX = x;
				scaleY = y;
				scaleZ = z;
			}

			void update(const liblas::Point &pt) {
				x = pt.GetX();
				y = pt.GetY();
				z = pt.GetZ();
				intensity = pt.GetIntensity();
				ret = pt.GetReturnNumber();
				numRets = pt.GetNumberOfReturns();
				cls = pt.GetClassification().GetClass();
				angle = pt.GetScanAngleRank();
			}

			void read(std::istream &str) {
				int32_t xx, yy, zz;
				str.read((char *) &xx, sizeof(int32_t));
				str.read((char *) &yy, sizeof(int32_t));
				str.read((char *) &zz, sizeof(int32_t));
				str.read((char *) &intensity, sizeof(uint16_t));
				str.read((char *) &ret, sizeof(uint16_t));
				str.read((char *) &numRets, sizeof(uint16_t));
				str.read((char *) &cls, sizeof(uint8_t));
				str.read((char *) &angle, sizeof(int8_t));
				x = (double) (xx * scaleX);
				y = (double) (yy * scaleY);
				z = (double) (zz * scaleZ);
			}

			void write(std::ostream &str) {
				int32_t xx = (int32_t) (x / scaleX);
				int32_t yy = (int32_t) (y / scaleY);
				int32_t zz = (int32_t) (z / scaleZ);
				str.write((char *) &xx, sizeof(int32_t));
				str.write((char *) &yy, sizeof(int32_t));
				str.write((char *) &zz, sizeof(int32_t));
				str.write((char *) &intensity, sizeof(uint16_t));
				str.write((char *) &ret, sizeof(uint16_t));
				str.write((char *) &numRets, sizeof(uint16_t));
				str.write((char *) &cls, sizeof(uint8_t));
				str.write((char *) &angle, sizeof(int8_t));
			}

			bool last() const {
				return numRets > 0 && ret == numRets;
			}

			bool first() const {
				return numRets > 0 && ret == 1;
			}

			bool intermediate() const {
				return numRets > 2 && ret > 1 && ret < numRets;
			}
			
			bool ground() const {
				return cls == 2;
			}

			bool single() const {
				return numRets == 1;
			}
		};

		class SortedPointStream {
		private:
			Bounds m_bounds;
			BoundsTracker m_boundsTracker;
			unsigned int m_pointCount;
			std::vector<unsigned int> m_rows;
			unsigned int m_row;

			void init(const std::list<std::string> &files, double blockSize);

		public:

			/**
			 * Constructs a new point stream on the given files.
			 *
			 * The block size gives the size, in map units, of a block,
			 * blocks are square and stored as individual files.
			 */
			SortedPointStream(const std::list<std::string> &files, double blockSize);

			~SortedPointStream();


			unsigned int pointCount() const;

			unsigned int rowCount() const;

			bool contains(double x, double y) const;

			bool contains(double x1, double y1, double x2, double y2) const;

			const Bounds& bounds() const;

			void next(std::list<LASPoint> &pts);
		};

	}
}


#endif