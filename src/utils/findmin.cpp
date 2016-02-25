#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Raster.hpp"
#include "Util.hpp"

int main(int argc, char **argv) {
	
	std::string infile = argv[1];

	try {
		Raster<float> inrast(infile);
		std::string proj;
		inrast.projection(proj);

		std::cerr << "Raster: " << inrast.cols() << ", " << inrast.rows() << std::endl;

		Grid<double> grid(inrast.cols(), inrast.rows());
		grid.fill(-9999);

		double nodata = inrast.nodata();
		for(int r = 0; r < inrast.rows(); ++r) {
			//std::cerr << "Row " << r << " of " << inrast.rows() << std::endl;
			for(int c = 0; c < inrast.cols(); ++c) {
				//std::cerr << c << ", " << r << std::endl;
				double z = inrast.getOrNodata(c, r);
				if(z == nodata) continue;
				bool min = true;
				for(int rr = r - 1; rr < r + 2 && !min; ++rr) {
					for(int cc = c - 1; cc < c + 2; ++cc) {
						double z0 = inrast.getOrNodata(cc, rr);
						if(z0 != nodata && z0 < z) {
							min = false;
							break;
						}
					}
				}

				if(min)
					grid(c, r, z);

			}
		}

		std::cout << std::setprecision(9) << "x,y,z" << std::endl;
		for(int r = 0; r < grid.rows(); ++r) {
			for(int c = 0; c < grid.cols(); ++c) {
				if(grid(c, r) > -9999)
					std::cout << inrast.toX(c) << "," << inrast.toY(r) << "," << grid(c, r) << std::endl;
			}
		}

	} catch(const char *e) {
		std::cerr << e << std::endl;
		return 1;
	}
	return 0;
}