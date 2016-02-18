#ifndef _INTERP_POINT_H_
#define _INTERP_POINT_H_

namespace interp {

	class InterpPoint {
	public:
		double x, y;
		double adj, base, change;
		InterpPoint() :
			x(nan("")),
			y(nan("")),
			adj(nan("")),
			base(nan("")),
			change(nan("")) {

		}
		InterpPoint(double x, double y) : InterpPoint() {
			this->x = x;
			this->y = y;
		}
		double diff() {
			return adj - base;
		}
		double resid() {
			return diff() - change;
		}
		bool equals(const InterpPoint &p) const {
			return x == p.x && y == p.y;
		}
		~InterpPoint() {
		}
	};

}

#endif
