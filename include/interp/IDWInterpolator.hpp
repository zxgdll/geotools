#ifndef _IDW_INTERPOLATOR_H_
#define _IDW_INTERPOLATOR_H_

#include "interp/Interpolator.hpp"

namespace interp {

	namespace idw {

		class IDWInterpolator : public interp::Interpolator {
		private:
			double i_exponent;
			unsigned int i_neighbours;
		public:

			IDWInterpolator(double exponent, unsigned int neighbours) {
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
