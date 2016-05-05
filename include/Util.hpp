#ifdef LIBLAS_HPP_INCLUDED
namespace las = liblas;
#endif

#include <set>
#include <list>
#include <sstream>
#include <algorithm>

#include "csv.h"


/**
 * Provides utility methods for working with LiDAR data.
 */
class Util {
public:

#ifdef LIBLAS_HPP_INCLUDED

	/**
	 * Compute the bounds of a LAS file using the Header and a double
	 * array to contain the result. If the dims argument is 2, the
	 * horizontal coords are considered and if it's 3, all coords are.
	 */
	static bool computeLasBounds(las::Header &hdr, std::vector<double> &bounds, int dims) {
		if((dims != 2 && dims != 3) || hdr.GetMaxX() == hdr.GetMinX() || hdr.GetMaxY() == hdr.GetMinY())
			return false;
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
		return true;
	}

	/**
	 * Compute the bounds of a LAS file using the Reader and a double
	 * array to contain the result. If the dims argument is 2, the
	 * horizontal coords are considered and if it's 3, all coords are.
	 * Use this method when the header bounds are bogus (it iterates over all the points.
	 */
	static bool computeLasBounds(las::Reader &rdr, std::vector<double> &bounds, int dims) {
		if(dims != 2 && dims != 3)
			return false;
		double x, y, z;
		while(rdr.ReadNextPoint()) {
			las::Point pt = rdr.GetPoint();
			x = pt.GetX();
			y = pt.GetY();
			if(x < bounds[0]) bounds[0] = x;
			if(y < bounds[1]) bounds[1] = y;
			if(x > bounds[2]) bounds[2] = x;
			if(y > bounds[3]) bounds[3] = y;
			if(dims == 3) {
				z = pt.GetZ();
				if(z < bounds[4]) bounds[4] = z;
				if(z > bounds[5]) bounds[5] = z;
			}
		}
		return true;
	}

#endif

	/**
	 * Expands the first bounds argument (an array of doubles) to contain
	 * the second (also an array of doubles). If the dims argument is 2,
	 * only horizontal dimensions are used. If 3, the z dimension is included.
	 */
	static void expand(std::vector<double> &bounds, std::vector<double> &extra, int dims) {
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
	static void snapBounds(std::vector<double> &bounds, double res, int dims) {
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
	static void printBounds(std::vector<double> &bounds, int dims) {
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
	static bool intersects(std::vector<double> &a, std::vector<double> &b, int dims) {
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
	static bool inBounds(double x, double y, std::vector<double> &bounds) {
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
	 * Split a comma-delimited string into a set of unique integers.
	 */
	static void intSplit(std::list<int> &values, const char *val) {
		std::stringstream ss(val);
		std::string item;
		while(std::getline(ss, item, ','))
			values.push_back(atoi(item.c_str()));
	}

	/**
	 * Split a comma-delimited string into a set of unique integers.
	 */
	static void intSplit(std::vector<int> &values, const char *str) {
		std::stringstream ss(str);
		std::string item;
		while(std::getline(ss, item, ','))
			values.push_back(atoi(item.c_str()));
	}

	/**
	 * Return true if the integer is in the set, or the set is empty.
	 */
	static bool inList(std::set<int> &values, int value) {
		return values.size() == 0 || values.find(value) != values.end();
	}

	static bool inList(std::vector<int> &values, int value) {
		return std::find(values.begin(), values.end(), value) != values.end();
	}

	/**
	 * Prints out a status message; a percentage representing current
	 * of total steps.
	 */
	static void status(int current, int total) {
		static int len = 0;
		float perc = (float) current / total * 100;
		if(len > 0) {
			for(int i = 0; i < len; ++i)
				std::cout << '\b';
		}
		std::stringstream ss;
		ss << perc << '%';
		std::cout << ss.str();
		len = ss.str().size();
	}

	static void copyfile(std::string &srcfile, std::string &dstfile) {
		std::ifstream src(srcfile.c_str(), std::ios::binary);
		std::ofstream dst(dstfile.c_str(), std::ios::binary);
		dst << src.rdbuf();
	}

	/**
	 * Load the samples from a csv file. The file must have x, y and z headers.
	 */
	static void loadXYZSamples(std::string &datafile, std::vector<std::tuple<double, double, double> > &samples) {
		io::CSVReader<3> in(datafile.c_str());
		in.read_header(io::ignore_extra_column, "x", "y", "z");
		double x, y, z;
		while(in.read_row(x, y, z))
			samples.push_back(std::make_tuple(x, y, z));
	}

	static void loadIDXYZSamples(std::string &datafile, std::vector<std::tuple<std::string, double, double, double> > &samples) {
		io::CSVReader<4> in(datafile.c_str());
		in.read_header(io::ignore_extra_column, "id", "x", "y", "z");
		std::string id;
		double x, y, z;
		while(in.read_row(id, x, y, z))
			samples.push_back(std::make_tuple(id, x, y, z));
	}


};
