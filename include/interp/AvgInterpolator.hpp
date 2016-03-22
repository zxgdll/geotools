#ifndef _AVG_INTERPOLATOR_
#define _AVG_INTERPOLATOR_

#include "interp/Interpolator.hpp"

namespace interp {

	namespace avg {

		class AvgInterpolator : public interp::Interpolator {
		private:
			unsigned int m_neighbours = 0;
		public:
			/**
			 * Create an AvgInterpolator with n neighbours. If n
			 * is <= 0, all data points are used and every cell has the
			 * same value. If n > 0, only the n nearest neighbours are
			 * averaged.
			 */
			AvgInterpolator(unsigned int neighbours) {
				m_neighbours = neighbours;
			}
			void interpolate(Raster<float> &out, std::list<InterpPoint> &samples);

		};

	}
}

#endif
