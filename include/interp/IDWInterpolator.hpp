#ifndef _IDW_INTERPOLATOR_H_
#define _IDW_INTERPOLATOR_H_

#include "interp/Interpolator.hpp"

namespace interp {

	namespace idw {

		class IDWInterpolator : public interp::Interpolator {
		private:
			double m_exponent;
			unsigned int m_neighbours;
		public:

			IDWInterpolator(double exponent, unsigned int neighbours) {
				if(exponent <= 0.0)
					throw "Please uses an exponent larger than zero.";
				m_exponent = exponent;
				m_neighbours = neighbours;
			}

			void interpolate(Raster<float> &out, std::list<InterpPoint> &samples);

		};

	}
}

#endif
