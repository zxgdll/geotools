#ifndef _IDW_INTERPOLATOR_H_
#define _IDW_INTERPOLATOR_H_

#include "interp/Interpolator.hpp"

namespace interp {

	namespace idw {

		class IDWInterpolator : public interp::Interpolator {
		private:
			double i_exponent;
		public:

			IDWInterpolator(double exponent) {
				i_exponent = exponent;
			}

			void interpolate(Raster<float> &out, std::list<InterpPoint> &samples);

		};

	}
}

#endif
