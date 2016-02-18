#ifndef _AVG_INTERPOLATOR_
#define _AVG_INTERPOLATOR_

#include "interp/Interpolator.hpp"

namespace interp {

	namespace avg {

		class AvgInterpolator : public interp::Interpolator {
		public:

			void interpolate(Raster<float> &out, std::list<InterpPoint> &samples);

		};

	}
}

#endif
