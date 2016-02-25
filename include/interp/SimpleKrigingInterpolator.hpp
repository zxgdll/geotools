#ifndef _SKRIGE_INTERPOLATOR_
#define _SKRIGE_INTERPOLATOR_

#include "interp/Interpolator.hpp"

namespace interp {

	namespace kriging {

		class VariogramPoint {
		private:
			double m_dist;
			double m_diff;
		public:
			VariogramPoint(double distance, double difference) {
				m_dist = distance;
				m_diff = difference;
			}
			double distance() {
				return m_dist;
			}
			double difference() {
				return m_diff;
			}
		};

		class SimpleKrigingInterpolator : public interp::Interpolator {
		private:
			int m_argc;
			char **m_argv;
			int argc() const {
				return m_argc;
			}
			char **argv() const {
				return m_argv;
			}
		public:

			SimpleKrigingInterpolator(int argc, char **argv) {
				m_argc = argc;
				m_argv = argv;
			}

			void computeVariogram(std::list<InterpPoint> &samples, std::list<VariogramPoint> &variogram);

			void showVariogram(std::list<VariogramPoint> &variogram);

			void interpolate(Raster<float> &out, std::list<InterpPoint> &samples);

		};

	}
}

#endif
