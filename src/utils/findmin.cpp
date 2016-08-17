#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Raster.hpp"
#include "util.hpp"

int g_min(int a, int b) {
	return a < b ? a : b;
}

int main(int argc, char **argv) {
	
	std::string infile = argv[1];

	try {
		Raster<float> inrast(infile);
		std::string proj;
		inrast.projection(proj);
		double nodata = inrast.nodata();

		std::cerr << "Raster: " << inrast.cols() << ", " << inrast.rows() << std::endl;

		Grid<double> grid(inrast.cols(), inrast.rows());
		grid.fill(nodata);

		int brows = 1024;
		float *block = (float *) malloc(sizeof(float) * inrast.cols() * brows);

		for(int br = 0; br < inrast.rows(); br += brows) {
			std::cerr << "Row " << br << " of " << inrast.rows() << std::endl;
			if(br + brows >= inrast.rows())
				brows = inrast.rows() - br;
			inrast.loadBlock(0, br, inrast.cols(), brows, block);
			for(int r = 1; r < brows - 1; ++r) {
				for(int c = 1; c < inrast.cols() - 1; ++c) {
					//std::cerr << c << ", " << r << std::endl;
					double z = block[r * inrast.cols() + c];
					if(z == nodata) continue;
					bool min = true;
					for(int rr = r - 1; min && (rr < r + 1); ++rr) {
						for(int cc = c - 1; min && (cc < c + 1); ++cc) {
							double z0 = block[rr * inrast.cols() + cc];
							if(z0 != nodata && z0 < z)
								min = false;
						}
					}

					if(min)
						grid(c, br + r, z);
				}
			}
		}

		std::cerr << "Writing..." << std::endl;
		std::cout << std::setprecision(9) << "x,y,z" << std::endl;
		for(int r = 0; r < grid.rows(); ++r) {
			for(int c = 0; c < grid.cols(); ++c) {
				if(grid(c, r) != nodata)
					std::cout << inrast.toX(c) << "," << inrast.toY(r) << "," << grid(c, r) << std::endl;
			}
		}

	} catch(const char *e) {
		std::cerr << e << std::endl;
		return 1;
	}
	return 0;
}