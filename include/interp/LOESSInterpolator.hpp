#ifndef _LOESS_INTERPOLATOR_H_
#define _LOESS_INTERPOLATOR_H_

#include "interp/Interpolator.hpp"

namespace interp {

    namespace loess {

        class LOESSInterpolator : public interp::Interpolator {
        private:
            double m_bandwidth;
            unsigned int m_degree;
        public:

            LOESSInterpolator(unsigned int degree, double bandwidth) {
                m_degree = degree;
                m_bandwidth = bandwidth;
            }

            void interpolate(Raster<float> &out, std::list<InterpPoint> &samples);

        };

    }
}

#endif
