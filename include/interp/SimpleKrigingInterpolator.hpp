#ifndef _SKRIGE_INTERPOLATOR_
#define _SKRIGE_INTERPOLATOR_

#include "interp/Interpolator.hpp"

namespace interp {

	namespace kriging {

		class VariogramPoint {
		private:
			double v_dist;
			double v_diff;
		public:
			VariogramPoint(double distance, double difference) {
				v_dist = distance;
				v_diff = difference;
			}
			double distance() {
				return v_dist;
			}
			double difference() {
				return v_diff;
			}
		};

		class SimpleKrigingInterpolator : public interp::Interpolator {
		private:
			int s_argc;
			char **s_argv;
			int argc() const {
				return s_argc;
			}
			char **argv() const {
				return s_argv;
			}
		public:

			SimpleKrigingInterpolator(int argc, char **argv) {
				s_argc = argc;
				s_argv = argv;
			}

			void computeVariogram(std::list<InterpPoint> &samples, std::list<VariogramPoint> &variogram);

			void showVariogram(std::list<VariogramPoint> &variogram);

			void interpolate(Raster<float> &out, std::list<InterpPoint> &samples);

		};

	}
}

#endif
