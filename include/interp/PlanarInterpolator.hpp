#ifndef _PLANAR_INTERPOLATOR_
#define _PLANAR_INTERPOLATOR_

#include "interp/Interpolator.hpp"

namespace interp {

    namespace planar {

        class PlanarInterpolator : public interp::Interpolator {
        public:

            void interpolate(Raster<float> &out, std::list<InterpPoint> &samples);

        };

    }
}

#endif
