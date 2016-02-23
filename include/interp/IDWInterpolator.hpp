#ifndef _IDW_INTERPOLATOR_H_
#define _IDW_INTERPOLATOR_H_

#include "interp/Interpolator.hpp"

namespace interp {

	namespace idw {

		class IDWInterpolator : public interp::Interpolator {
		private:
			double i_exponent;
			int i_neighbours;
		public:

			IDWInterpolator(double exponent, int neighbours) {
				if(neighbours < 0)
					throw "Invalid neighbours count. Must be >= 0.";
				if(exponent <= 0.0)
					throw "Please uses an exponent larger than zero.";
				i_exponent = exponent;
				i_neighbours = neighbours;
			}

			void interpolate(Raster<float> &out, std::list<InterpPoint> &samples);

		};

	}
}

#endif
