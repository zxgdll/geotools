#ifndef __SIMPLEGEOM_HPP__
#define __SIMPLEGEOM_HPP__

#include "geos/geom/GeometryFactory.h"
#include "geos/geom/Polygon.h"
#include "geos/geom/Geometry.h"
#include "geos/util/GeometricShapeFactory.h"
#include "geos/geom/Point.h"
#include "geos/geom/LinearRing.h"
#include "geos/geom/CoordinateSequence.h"
#include "geos/geom/CoordinateSequenceFactory.h"

namespace util = geotools::util;

namespace geotools {
	namespace geom {

		const geos::geom::GeometryFactory gf;
		geos::util::GeometricShapeFactory gsf(&gf);
		
		class SimpleGeom {
		private:

		public:

			static std::unique_ptr<geos::geom::Point> createPoint(double x, double y, double z) {
				const geos::geom::Coordinate c(x, y, z);
				geos::geom::Point *p = gf.createPoint(c);
				std::unique_ptr<geos::geom::Point> pt(p);
				return std::move(pt);
			}

			static std::unique_ptr<geos::geom::Polygon> createCircle(const geos::geom::Coordinate &coord, double radius) {
				gsf.setBase(coord);
				gsf.setSize(radius * 2.0);
				std::unique_ptr<geos::geom::Polygon> circle(gsf.createCircle());
				return std::move(circle);
			}

			static std::unique_ptr<geos::geom::Polygon> createPolygon(const std::list<geos::geom::Coordinate> &coords, 
					const std::list<std::list<geos::geom::Coordinate> > &holes) {
				geos::geom::CoordinateSequence *cs = gf.getCoordinateSequenceFactory()->create(nullptr, 3);
				for(const geos::geom::Coordinate &coord : coords)
					cs->add(coord, false);
				// TODO: Holes
				geos::geom::LinearRing ring(cs, &gf);
				geos::geom::Polygon *ply = gf.createPolygon(&ring, nullptr);
				std::unique_ptr<geos::geom::Polygon> poly(ply);
				return std::move(poly);
			}

			static std::unique_ptr<geos::geom::Polygon> createPolygon(const std::list<geos::geom::Coordinate> &coords) {
				std::list<std::list<geos::geom::Coordinate> > holes;
				return std::move(createPolygon(coords, holes));
			}

			static std::unique_ptr<geos::geom::Polygon> createPolygon(const geotools::util::Bounds &bounds) {
				std::list<geos::geom::Coordinate> coords {
					geos::geom::Coordinate(bounds.minx(), bounds.miny()),
					geos::geom::Coordinate(bounds.maxx(), bounds.miny()),
					geos::geom::Coordinate(bounds.maxx(), bounds.maxy()),
					geos::geom::Coordinate(bounds.minx(), bounds.maxy()),
					geos::geom::Coordinate(bounds.minx(), bounds.miny())
				};
				return std::move(createPolygon(coords));
			}

		};

	}
}

#endif