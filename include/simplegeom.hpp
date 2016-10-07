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

using namespace geos::geom;
using namespace geos::util;

namespace geotools {
	namespace geom {

		static const GeometryFactory *gf = GeometryFactory::getDefaultInstance();
		static GeometricShapeFactory *gsf = new GeometricShapeFactory(gf);

		class SimpleGeom {
		public:

			static std::unique_ptr<geos::geom::Point> createPoint(double x, double y, double z) {
				const Coordinate c(x, y, z);
				geos::geom::Point *p = gf->createPoint(c);
				std::unique_ptr<geos::geom::Point> pt(p);
				return std::move(pt);
			}

			static std::unique_ptr<Polygon> createCircle(const Coordinate &coord, double radius) {
				gsf->setBase(coord);
				gsf->setSize(radius * 2.0);
				std::unique_ptr<Polygon> circle(gsf->createCircle());
				return std::move(circle);
			}

			static std::unique_ptr<Polygon> createPolygon(const std::list<Coordinate> &coords, 
					const std::list<std::list<Coordinate> > &holes) {
				CoordinateSequence *cs = gf->getCoordinateSequenceFactory()->create(nullptr, 3);
				for(const Coordinate &coord : coords)
					cs->add(coord, false);
				cs->add(coords.front(), false);
				// TODO: Holes
				LinearRing *ring = new LinearRing(cs, gf);
				Polygon *ply = gf->createPolygon(ring, nullptr);
				std::unique_ptr<Polygon> poly(ply);
				return std::move(poly);
			}

			static std::unique_ptr<Polygon> createPolygon(const std::list<Coordinate> &coords) {
				std::list<std::list<Coordinate> > holes;
				return std::move(createPolygon(coords, holes));
			}

			static std::unique_ptr<Polygon> createPolygon(const geotools::util::Bounds &bounds) {
				std::list<Coordinate> coords {
					Coordinate(bounds.minx(), bounds.miny()),
					Coordinate(bounds.maxx(), bounds.miny()),
					Coordinate(bounds.maxx(), bounds.maxy()),
					Coordinate(bounds.minx(), bounds.maxy()),
					Coordinate(bounds.minx(), bounds.miny())
				};
				return std::move(createPolygon(coords));
			}

		};

	}
}

#endif