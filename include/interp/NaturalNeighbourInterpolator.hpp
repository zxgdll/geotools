#ifndef _NATN_INTERPOLATOR_
#define _NATN_INTERPOLATOR_

#include "interp/Interpolator.hpp"

namespace interp {

	namespace naturalneighbour {

		class NaturalNeighbourInterpolator : public interp::Interpolator {
		public:

			void interpolate(Raster<float> &out, std::list<InterpPoint> &samples);

		};

	}
}

#endif
