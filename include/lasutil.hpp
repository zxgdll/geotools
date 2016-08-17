#ifndef __LASUTIL_HPP__
#define __LASUTIL_HPP__

#include <liblas/liblas.hpp>

#include "util.hpp"

namespace geotools {

	namespace util {

		class LasUtil {
		public:

                        /**
                         * Compute the bounds of a LAS file using the Header and a double
                         * array to contain the result. If the dims argument is 2, the
                         * horizontal coords are considered and if it's 3, all coords are.
                         */
                        static bool computeLasBounds(liblas::Header &hdr, geotools::util::Bounds &bounds, int dims = 2);

                        /**
                         * Compute the bounds of a LAS file using the Reader and a double
                         * array to contain the result. If the dims argument is 2, the
                         * horizontal coords are considered and if it's 3, all coords are.
                         * Use this method when the header bounds are bogus (it iterates over all the points.
                         */
                        static bool computeLasBounds(liblas::Reader &rdr, geotools::util::Bounds &bounds, int dims = 2);


		};

	} // util

} // geotools

#endif
