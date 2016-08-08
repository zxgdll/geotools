#include "geotools.h"

#include "trees.hpp"

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
			} else if(arg == "-v") {
				_loglevel(LOG_TRACE);
			}
		}

		trees::treetops(inraster, topshp, window);

	} catch(const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;
}


