#include <string>
#include <vector>
#include <map>

#include "omp.h"

#include "geotools.h"
#include "treetops.hpp"

void usage() {
	std::cerr << "Usage: treecrowns <options>\n"
			<< "This program finds tree tops by locating the maximum value in a \n"
			<< "window as it moves across a raster. The tops are just the maxima at \n"
			<< "the center of the window.\n"
			<< " -i  <filename>     The input raster; a LiDAR-derived canopy height model.\n"
			<< " -t  <filename>     The treetop vector file. An sqlite file.\n"
			<< " -w  <int>          Window size. Will be bumped up to the next odd value if \n"
			<< "                    even given. Min 3.\n"
			<< " -sf <filename>     If the CHM is to be smoothed, enter the filename of the\n"
			<< "                    smoothed raster.\n"
			<< " -ss <float>        The std. deviation value for gaussian smoothing (if -sf is\n"
			<< " -sw <int>          The smoothing window size. 3 or greater; odd.\n"
			<< "                    provided). Default 0.8408964\n"
			<< " -cr <filename>     The crowns raster. If this is not provided, crowns are not\n" 
			<< "                    produced.\n"
			<< " -cv <filename>     The crowns vector. A vectorized version of the crown \n"
			<< "                    raster. Optional.\n"
			<< " -cm <float>        The minimum height of pixels to consider for inclusion \n"
			<< "                    in a crown. Default 4.\n"
			<< " -ct <float>        The minimum height of pixels, as a proportion of the \n"
			<< "                    treetop height\n"
			<< "                    to consider for inclusion in a crown. 0 < n < 1; default 0.65.\n"
			<< " -cd <float>        The radius beyond which pixels will not be considered for \n"
			<< "                    inclusion in\n"
			<< "                    in a crown. This value is also used for spatially partitioning \n"
			<< "                    the delineation work, so a reasonable value should be selected. \n"
			<< "                    Default 5; max 100m.\n"
			<< " -tm <float>        The minimum height of pixels to consider for selection as a \n"
			<< "                    tree top. Default 4.\n"
			<< " -d8                Use D8 search for delineating crowns, rather than D4, which \n"
			<< "                    is the default.\n"
			<< " -threads           The number of threads to used for processing. Default 1.\n"
			<< " -s                 The SRID of the treetops output dataset.\n";
			
}

int main(int argc, char **argv) {

	using namespace geotools::trees;
	using namespace geotools::trees::util;
	using namespace geotools::trees::config;

	try {

		TreeTopConfig ttConfig;
		std::string inraster;  // Input raster.
		std::string topsvect;  // Treetops output.
		std::string crownrast;
		std::string crownvect;
		std::string smoothed;
		double threshold = 0.65;
		double radius = 5;
		double cminHeight = 4;
		bool d8 = false;
		int threads = 1;

		int i = 1;
		for(; i < argc; ++i) {
			std::string arg = argv[i];
			if(arg == "-i") {
				ttConfig.setInputFilename(argv[++i]);
			} else if(arg == "-t") {
				ttConfig.setOutputFilename(argv[++i]);
			} else if(arg == "-s") {
				ttConfig.setSRID(atoi(argv[++i]));
			} else if(arg == "-cr") {
				crownrast = argv[++i];
			} else if(arg == "-cv") {
				crownvect = argv[++i];
			} else if(arg == "-cm") {
				cminHeight = atof(argv[++i]);
			} else if(arg == "-ct") {
				threshold = atof(argv[++i]);
			} else if(arg == "-cd") {
				radius = atof(argv[++i]);
			} else if(arg == "-tm") {
				ttConfig.setMinHeight(atof(argv[++i]));
			} else if(arg == "-w") {
				ttConfig.setSearchWindow(atoi(argv[++i]));
			} else if(arg == "-v") {
				g_loglevel(G_LOG_TRACE);
			} else if(arg == "-sf") {
				ttConfig.setSmoothedFilename(argv[++i]);
			} else if(arg == "-ss") {
				ttConfig.setSmoothingSigma(atof(argv[++i]));
			} else if(arg == "-sw") {
				ttConfig.setSmoothingWindow(atoi(argv[++i]));
			} else if(arg == "-d8") {
				d8 = true;
			} else if(arg == "-threads") {
				threads = atoi(argv[++i]);
			}
		}

		if(threads <= 0)
			g_argerr("Invalid number of threads: " << threads);
		omp_set_dynamic(0);
		omp_set_num_threads(threads);

		bool crowns = true;
		if(crownrast.empty()) {
			g_warn("A crown raster file is required for producing crowns. Crowns will not be computed.");
			crowns = false;
		}

		Trees tu;

		// Create tree tops.
		tu.treetops(ttConfig);

		// Substitute the input raster for the smoothed raster (because there isn't one)
		if(ttConfig.smoothedFilename.empty())
			ttConfig.setSmoothedFilename(ttConfig.inputFilename);

		// Create crowns if desired.
		if(crowns)
			tu.treecrowns(ttConfig.smoothedFilename, ttConfig.outputFilename, crownrast, crownvect, threshold, radius, cminHeight, d8);

	} catch(const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;
}


