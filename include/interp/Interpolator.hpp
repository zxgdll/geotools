#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include <list>
#include "Raster.hpp"
#include "interp/InterpPoint.hpp"

namespace interp {

    class Interpolator {
    public:

        /**
         * Interpolates the samples into the given raster.
         */
        virtual void interpolate(Raster<float> &out, std::list<InterpPoint> &samples) = 0;

        virtual ~Interpolator() {
        }

    };

}

#endif
