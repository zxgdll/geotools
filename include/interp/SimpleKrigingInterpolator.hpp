#ifndef _SKRIGE_INTERPOLATOR_
#define _SKRIGE_INTERPOLATOR_

#include "interp/Interpolator.hpp"

namespace interp {

	namespace kriging {

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

			void showVariogram(std::list<InterpPoint> &samples);

			void interpolate(Raster<float> &out, std::list<InterpPoint> &samples);

		};

	}
}

#endif
