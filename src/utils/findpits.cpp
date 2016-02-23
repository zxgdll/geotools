#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "Raster.hpp"
#include "Util.hpp"

int _max(int a, int b) {
	return a > b ? a : b;
}

int _min(int a, int b) {
	return a < b ? a : b;
}

int main(int argc, char **argv) {
	
	std::string infile = argv[1];
	double resolution = atoi(argv[2]);

	try {
		Raster<float> inrast(infile);
		std::string proj;
		inrast.projection(proj);

		int block = (int) (resolution / inrast.resolution());
		std::cerr << "Block size: " << block << std::endl;

		Grid<double> grid(inrast.cols() / block, inrast.rows() / block);
		grid.fill(-9999);

		for(int r = 0; r < inrast.rows(); r += block) {
			std::cerr << "Row " << r << " of " << inrast.rows() << std::endl;
			for(int c = 0; c < inrast.cols(); c += block) {
				if(!inrast.isValid(c, r)) continue;
				double z = inrast.get(c, r);
				int lc, lr;
				for(int rr = r; rr < _min(r + block, inrast.rows()); ++rr) {
					for(int cc = c; cc < _min(c + block, inrast.cols()); ++cc) {
						if(cc == c || rr == r || !inrast.isValid(cc, rr)) continue;
						double z0 = inrast.get(cc, rr);
						if(z0 < z) {
							z = z0;
							lc = cc;
							lr = rr;
						}
					}
				}
				grid(lc / block, lr / block, z);
			}
		}

		std::cout << std::setprecision(9);
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