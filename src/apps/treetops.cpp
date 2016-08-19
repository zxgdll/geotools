#include <string>
#include <vector>
#include <map>

#include "omp.h"

#include "geotools.h"
#include "trees.hpp"

void usage() {
	std::cerr << "Usage: treecrowns <options>\n"
			<< "This program finds tree tops by locating the maximum value in a \n"
			<< "window as it moves across a raster. The tops are just the maxima at \n"
			<< "the center of the window.\n"
			<< " -i  <filename>     The input raster; a LiDAR-derived canopy height model.\n"
			<< " -t  <filename>     The treetop vector file. An sqlite file.\n"
			<< " -w  <size>         Window size. Will be bumped up to the next odd value if even given.\n"
			<< " -sf <filename>     If the CHM is to be smoothed, enter the filename of the smoothed raster.\n";
			
}

int main(int argc, char **argv) {

	using namespace geotools::trees;
	using namespace geotools::trees::util;

	try {
		std::string inraster;  // Input raster.
		std::string topsvect;  // Treetops output.
		std::string crownrast;
		std::string crownvect;
		std::string smoothed;
		int window = 3;
		double tminHeight = 4;
		double threshold = 0.65;
		double radius = 5;
		double cminHeight = 4;

		int i = 1;
		for(; i < argc; ++i) {
			std::string arg = argv[i];
			if(arg == "-i") {
				inraster = argv[++i];
			} else if(arg == "-t") {
				topsvect = argv[++i];
			} else if(arg == "-cr") {
				crownrast = argv[++i];
			} else if(arg == "-cv") {
				crownvect = argv[++i];
			} else if(arg == "-w") {
				window = atoi(argv[++i]);
			} else if(arg == "-v") {
				g_loglevel(G_LOG_TRACE);
			} else if(arg == "-sf") {
				smoothed.assign(argv[++i]);
			} else if(arg == "-threads") {
				int t = atoi(argv[++i]);
				if(t <= 0)
					g_argerr("Invalid number of threads: " << t);
				omp_set_dynamic(0);
				omp_set_num_threads(t);
			}
		}

		bool crowns = true;
		if(crownrast.empty()) {
			g_warn("A crown raster file is required for producing crowns. Crowns will not be computed.");
			crowns = false;
		}

		TreeUtil tu;

		// Create a smoothed raster.
		if(!smoothed.empty()) {
			tu.smooth(inraster, smoothed);
			inraster.assign(smoothed);
		}

		// Create tree tops.
		std::vector<std::unique_ptr<Top> > tops;
		tu.treetops(inraster, topsvect, window, tminHeight, crowns ? &tops : nullptr);

		// Create crowns if desired.
		if(crowns)
			tu.treecrowns(inraster, tops, crownrast, crownvect, threshold, radius, cminHeight);

void treecrowns(const std::string &inraster, const std::vector<std::unique_ptr<geotools::trees::util::Top> > &tops, 
                                const std::string &crownsrast, const std::string &crownsvect, double threshold, double radius, 
                                double minHeight, bool d8 = false);

	} catch(const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;
}


