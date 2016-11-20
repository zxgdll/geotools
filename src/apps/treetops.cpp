#include <string>
#include <vector>
#include <map>

#include "omp.h"

#include "geotools.h"
#include "treetops.hpp"

#ifdef WITH_GUI
#include "treetops_ui.hpp"
#endif

void usage() {
	std::cerr << "Usage: treecrowns <options>\n"
			<< "This program finds tree tops by locating the maximum value in a \n"
			<< "window as it moves across a raster. The tops are just the maxima at \n"
			<< "the center of the window.\n"
			<< " -si <filename>     The input raster to be smoothed.\n"
			<< " -ss <int>          The standard deviation for Gaussian smoothing.\n"
			<< " -sw <int>          The window size for Gaussian smoothing. Will be \n"
			<< "                    bumped up to the next odd value if even given. Minimum 3.\n"
			<< " -so <filename>     The smoothed output raster.\n\n"

			<< " -ti <filename>     The input raster for treetop detection; a (smoothed)\n"
			<< "                    LiDAR-derived canopy height model.\n"
			<< " -tm <float>        The minimum height of pixels to consider for selection as a \n"
			<< "                    tree top. Default 4.\n"
			<< " -tw <int>          Window size. Will be bumped up to the next odd value if \n"
			<< "                    even given. Minimum 3.\n\n"
			<< " -to <filename>     The treetop vector file. An sqlite file.\n\n"

			<< " -ci <filename>     The input raster for crown delineation; a (smoothed)\n"
			<< "                    LiDAR-derived canopy height model.\n"
			<< " -ct <filename>     The treetop vector file. An sqlite file.\n\n"
			<< " -cm <float>        The minimum height of pixels to consider for inclusion \n"
			<< "                    in a crown. Default 4.\n"
			<< " -cf <float>        The minimum height of pixels, as a proportion of the \n"
			<< "                    treetop height to consider for inclusion in a crown.\n"
			<< "                    0 < n < 1; default 0.65.\n"
			<< " -cd <float>        The radius beyond which pixels will not be considered for \n"
			<< "                    inclusion in a crown.\n"
			<< " -cr <filename>     The crowns raster. Each crown is represented with pixels\n"
			<< "                    having a unique ID.\n"
			<< " -cv <filename>     The crowns vector. A vectorized version of the crown \n"
			<< "                    raster. Optional.\n\n"

			<< " -threads           The number of threads to used for processing. Default 1.\n"			
			<< " -v                 Print debug messages.\n";
}

int runWithGui(int argc, char **argv) {
#ifdef WITH_GUI
	QApplication q(argc, argv);
	QWidget *w = new QWidget();
	geotools::ui::TreetopsForm f;
	f.setupUi(w);
	w->show();
	return q.exec();
#else
	std::cerr << "GUI not enabled." << std::endl;
	return 1;
#endif
}

int main(int argc, char **argv) {

	using namespace geotools::treetops;
	using namespace geotools::treetops::util;
	using namespace geotools::treetops::config;

	try {

		TreetopsConfig config;

		int i = 1;
		for(; i < argc; ++i) {
			std::string arg = argv[i];

			if(arg == "-si") {
				config.smoothInputFile = argv[++i];
				config.doSmoothing = true;
			} else if(arg == "-ss") {
				config.smoothStdDev = atof(argv[++i]);
				config.doSmoothing = true;
			} else if(arg == "-so") {
				config.smoothOutputFile = argv[++i];
				config.doSmoothing = true;

			} else if(arg == "-ti") {
				config.topsInputFile = argv[++i];
				config.doTops = true;
			} else if(arg == "-tm") {
				config.topsMinHeight = atof(argv[++i]);
				config.doTops = true;
			} else if(arg == "-tw") {
				config.topsWindowSize = atoi(argv[++i]);
				config.doTops = true;
			} else if(arg == "-to") {
				config.topsOutputFile = argv[++i];
				config.doTops = true;

			} else if(arg == "-ci") {
				config.crownsInputFile = argv[++i];
				config.doCrowns = true;
			} else if(arg == "-ct") {
				config.crownsTreetopsFile = argv[++i];
				config.doCrowns = true;
			} else if(arg == "-cm") {
				config.crownsMinHeight = atof(argv[++i]);
				config.doCrowns = true;
			} else if(arg == "-cf") {
				config.crownsHeightFraction = atof(argv[++i]);
				config.doCrowns = true;
			} else if(arg == "-cd") {
				config.crownsRadius = atof(argv[++i]);
				config.doCrowns = true;
			} else if(arg == "-cr") {
				config.crownsOutputRaster = argv[++i];
				config.doCrowns = true;
			} else if(arg == "-cv") {
				config.crownsOutputVector = argv[++i];
				config.doCrowns = true;
			} else if(arg == "-threads") {
				config.threads = atoi(argv[++i]);
				config.doCrowns = true;
			}
		}

		Treetops tu;

		// Create tree tops.
		if(config.doSmoothing)
			tu.smooth(config);

		// Create treetops.
		if(config.doTops)
			tu.treetops(config);

		// Create crowns.
		if(config.doCrowns)
			tu.treecrowns(config);

		if(!config.doTops && !config.doSmoothing && !config.doCrowns)
			return runWithGui(argc, argv);

	} catch(const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;
}


