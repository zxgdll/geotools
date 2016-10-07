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

			static geos::geom::Point* createPoint(double x, double y, double z) {
				const Coordinate c(x, y, z);
				return gf->createPoint(c);
			}

			static Polygon* createCircle(const Coordinate &coord, double radius) {
				gsf->setBase(coord);
				gsf->setSize(radius * 2.0);
				return gsf->createCircle();
			}

			static Polygon* createPolygon(const std::list<Coordinate> &coords, 
					const std::list<std::list<Coordinate> > &holes) {
				CoordinateSequence *cs = gf->getCoordinateSequenceFactory()->create(nullptr, 3);
				for(const Coordinate &coord : coords)
					cs->add(coord, false);
				cs->add(coords.front(), false);
				// TODO: Holes
				LinearRing *ring = new LinearRing(cs, gf);
				return gf->createPolygon(ring, nullptr);
			}

			static Polygon* createPolygon(const std::list<Coordinate> &coords) {
				std::list<std::list<Coordinate> > holes;
				return createPolygon(coords, holes);
			}

			static Polygon* createPolygon(const geotools::util::Bounds &bounds) {
				std::list<Coordinate> coords {
					Coordinate(bounds.minx(), bounds.miny()),
					Coordinate(bounds.maxx(), bounds.miny()),
					Coordinate(bounds.maxx(), bounds.maxy()),
					Coordinate(bounds.minx(), bounds.maxy()),
					Coordinate(bounds.minx(), bounds.miny())
				};
				return createPolygon(coords);
			}

		};

	}
}

#endif