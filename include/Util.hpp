#include <liblas/liblas.hpp>

namespace las = liblas;

/**
 * Provides utility methods for working with LiDAR data.
 */
class Util {
public:

	/**
	 * Compute the bounds of a LAS file using the Header and a double
	 * array to contain the result. If the dims argument is 2, the
	 * horizontal coords are considered and if it's 3, all coords are.
	 */
	static void computeLasBounds(las::Header &hdr, double *bounds, int dims) {
		if(dims == 2 or dims == 3) {
			if(hdr.GetMinX() < bounds[0]) bounds[0] = hdr.GetMinX();
			if(hdr.GetMinY() < bounds[1]) bounds[1] = hdr.GetMinY();
			if(hdr.GetMaxX() > bounds[2]) bounds[2] = hdr.GetMaxX();
			if(hdr.GetMaxY() > bounds[3]) bounds[3] = hdr.GetMaxY();
		}
		if(dims == 3) {
			if(hdr.GetMinZ() < bounds[4]) bounds[5] = hdr.GetMinZ();
			if(hdr.GetMaxZ() > bounds[5]) bounds[5] = hdr.GetMaxZ();
		}
	}

	/**
	 * Expands the first bounds argument (an array of doubles) to contain
	 * the second (also an array of doubles). If the dims argument is 2,
	 * only horizontal dimensions are used. If 3, the z dimension is included.
	 */
	static void expand(double *bounds, double *extra, int dims) {
		if(dims == 2 || dims == 3) {
			if(extra[0] < bounds[0]) bounds[0] = extra[0];
			if(extra[1] < bounds[1]) bounds[1] = extra[1];
			if(extra[2] > bounds[2]) bounds[2] = extra[2];
			if(extra[3] > bounds[3]) bounds[3] = extra[3];
		}
		if(dims == 3) {
			if(extra[4] < bounds[4]) bounds[4] = extra[4];
			if(extra[5] > bounds[5]) bounds[5] = extra[5];
		}
	}

	/**
	 * Snaps the bounds to the given resolution. If the dims
	 * argument is 2, considers x and y. If 3, considers z.
	 */
	static void snapBounds(double *bounds, double res, int dims) {
		if(dims == 2 || dims == 3) {
			bounds[0] = floor(bounds[0] / res) * res;
			bounds[1] = floor(bounds[1] / res) * res;
			bounds[2] = ceil(bounds[2] / res) * res;
			bounds[3] = ceil(bounds[3] / res) * res;
		}
		if(dims == 3) {
			bounds[4] = floor(bounds[4] / res) * res;
			bounds[5] = ceil(bounds[5] / res) * res;
		}
	}

	/**
	 * Prints the bounds to stderr. If dims is 2, prints x and y.
	 * If 3 prints z.
	 */
	static void printBounds(double *bounds, int dims) {
		if(dims == 2 || dims == 3) {
			std::cerr << "Bounds: " << bounds[0] << "," << bounds[1] << "," << bounds[2] << "," << bounds[3];
			if(dims == 3)
				std::cerr << "," << bounds[4] << "," << bounds[5];
			std::cerr << std::endl;
		}
	}

	/**
	 * Returns true if the two bounds intersect. The dims argument
	 * determines whether the comparison is 2d (2) or 3d (3).
	 */
	static bool intersects(double *a, double *b, int dims) {
		if(dims == 2)
			return !(a[0] > b[2] || a[1] > b[3] || a[2] < b[0] || a[3] < b[1]);
		if(dims == 3)
			return !(a[0] > b[2] || a[1] > b[3] || a[2] < b[0] || a[3] < b[1] || a[4] > b[5] || a[5] < b[4]);
		return false;
	}

	/**
	* Returns true if the given coordinate is inside the bounds.
	* Coordinates are assumed to have the same CRS as the bounds.
	*/
	static bool inBounds(double x, double y, double *bounds) {
		return x >= bounds[0] && x < bounds[2] && y >= bounds[1] && y < bounds[3];
	}

	/**
	* Returns true if the given coordinate is inside the bounds.
	* Coordinates are assumed to have the same CRS as the bounds.
	*/
	static bool inBounds(double x, double y, double z, double *bounds) {
		return x >= bounds[0] && x < bounds[2] && y >= bounds[1] && y < bounds[3] && z >= bounds[4] && z < bounds[5];
	}

	/**
	 * Split a comma-delimited string into a set of unique integers.
	 */
	static void intSplit(std::set<int> &values, const char *str) {
		std::stringstream ss(str);
		std::string item;
		while(std::getline(ss, item, ','))
			values.insert(atoi(item.c_str()));
	}

	/**
	 * Return true if the integer is in the set, or the set is empty.
	 */
	static bool inList(std::set<int> &values, int value) {
		return values.size() == 0 || values.find(value) != values.end();
	}

};
