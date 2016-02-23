#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Raster.hpp"
#include "Util.hpp"

int main(int argc, char **argv) {
	
	std::string infile = argv[1];
	double resolution = atoi(argv[2]);

	try {
		Raster<float> inrast(infile);
		std::string proj;
		inrast.projection(proj);

		std::cerr << "Raster: " << inrast.cols() << ", " << inrast.rows() << std::endl;

		int block = (int) (resolution / inrast.resolution());
		std::cerr << "Block size: " << block << std::endl;

		Grid<double> grid(inrast.cols() / block + 1, inrast.rows() / block + 1);
		grid.fill(-9999);

		for(int r = 0; r < inrast.rows(); r += block) {
			//std::cerr << "Row " << r << " of " << inrast.rows() << std::endl;
			for(int c = 0; c < inrast.cols(); c += block) {
				//std::cerr << c << ", " << r << std::endl;
				if(!inrast.isValid(c, r)) continue;
				double z = inrast.get(c, r);
				int lc, lr;
				bool cont = true;
				do {
					cont = false;
					for(int rr = r - 1; rr < r + 2; ++rr) {
						for(int cc = c - 1; cc < c + 2; ++cc) {
							if(cc < 0 || rr < 0 || cc == c || rr == r || rr >= inrast.rows() || cc >= inrast.cols() || !inrast.isValid(cc, rr)) continue;
							double z0 = inrast.get(cc, rr);
							if(z0 < z) {
								z = z0;
								lc = cc;
								lr = rr;
								cont = true;
								//std::cerr << " -> " << lc << ", " << lr << ", " << z << std::endl;
							}
						}
					}
				} while(cont);

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