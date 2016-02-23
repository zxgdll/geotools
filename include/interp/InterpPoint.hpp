#ifndef _INTERP_POINT_H_
#define _INTERP_POINT_H_

namespace interp {

	class InterpPoint {
	public:
		double x, y, z, resid;
		InterpPoint() :
			x(nan("")),
			y(nan("")),
			z(nan("")),
			resid(nan("")) {

		}
		InterpPoint(double x, double y) : InterpPoint() {
			this->x = x;
			this->y = y;
		}
		InterpPoint(double x, double y, double z) : InterpPoint(x, y) {
			this->z = z;
		}
		bool equals(const InterpPoint &p) const {
			return x == p.x && y == p.y;
		}
		~InterpPoint() {
		}
	};

}

#endif
