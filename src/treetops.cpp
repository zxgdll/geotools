/*
 * treecrowns.cpp
 *
 *  Created on: May 3, 2016
 *      Author: rob
 */

#include <queue>
#include <iostream>

#include "Raster.hpp"
#include "Vector.hpp"

class Top {
public:
	size_t m_id;
	double m_x;
	double m_y;
	double m_z;
	Top(size_t id, double x, double y, double z) :
		m_id(id),
		m_x(x), m_y(y), m_z(z) {
	}
	Top(const Top &t) :
		Top(t.m_id, t.m_x, t.m_y, t.m_z) {
	}
	Top() :
		Top(0, 0, 0, 0) {
	}
};

bool findTop(MemRaster<float> &raster, int window) {
	float max = raster.get(window / 2, window / 2);
	if(max == raster.nodata())
		return false;
	for(int r = 0; r < raster.rows(); ++r) {
		for(int c = 0; c < raster.cols(); ++c) {
			if(c != window / 2 && r != window / 2 && raster.get(c, r) >= max)
				return false;
		}
	}
	return true;
}

/**
*/
 void treetops(std::string &inraster, std::string &topshp, int window) {

	if(inraster.empty())
		_argerr("Input raster cannot be empty.");
	if(topshp.empty())
		_argerr("The treetop output name cannot be empty.");
	if(window < 3)
		_argerr("A window size less than 3 makes no sense.");
	if(window % 2 == 0) {
		_warn("Window is " << window << ". Bumping up to " << (window + 1));
		window++;
	}

	// Input raster.
	Raster<float> raster(inraster);
	MemRaster<float> blk(window, window);
	blk.nodata(raster.nodata());

	std::map<size_t, Top*> tops;

	raster.setCacheSize(window);

	// Run along columns because that's usually the way gdal blocks are oriented.
	for(int c = 0; c < raster.cols() - window - 1; ++c) {
		for(int r = 0; r < raster.rows() - window - 1; ++r) {
			raster.readBlock(c, r, blk);
			if(findTop(blk, window)) {
				size_t id = ((size_t) c << 32) | r;
				if(tops.find(id) == tops.end()) {
					int off = window / 2;
					tops[id] = new Top(id, 
						raster.toX(c + off) + raster.resolutionX() / 2.0, 
						raster.toY(r + off) + raster.resolutionY() / 2.0, 
						blk.get(off, off)
					);
				}
			}
		}
	}

	std::stringstream out;
	out << std::setprecision(12);
	out << "id,x,y,z" << std::endl;
	for(auto it = tops.begin(); it != tops.end(); ++it) {
		Top *t = it->second;
		out << t->m_id << "," << t->m_x << "," << t->m_y << "," << t->m_z << std::endl;
		delete t;
	}
	std::cout << out.str();
}

void usage() {
	std::cerr << "Usage: treecrowns <options>\n"
			<< "This program finds tree tops by locating the maximum value in a \n"
			<< "window as it moves across a raster. The tops are just the maxima at \n"
			<< "the center of the window.\n"
			<< " -i -- The input raster; a LiDAR-derived canopy height model.\n"
			<< " -t -- The treetop vector file. A shapefile.\n"
			<< " -w -- Window size. Will be bumped up to the next odd value if even given.\n";
}

int main(int argc, char **argv) {

	try {
		std::string inraster; // Input raster.
		std::string topshp;   // Treetops output.
		int window = 0;

		int i = 1;
		for(; i < argc; ++i) {
			std::string arg = argv[i];
			if(arg == "-i") {
				inraster = argv[++i];
			} else if(arg == "-t") {
				topshp = argv[++i];
			} else if(arg == "-w") {
				window = atoi(argv[++i]);
			}
		}

		treetops(inraster, topshp, window);

	} catch(const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;
}



