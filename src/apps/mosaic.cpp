/**
 * "Feathers" the edges of a data region (i.e. not nodata) to the specified
 * distance in map units, using the specified curve.
 */
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <omp.h>

#include "geotools.h"
#include "mosaic.hpp"

#ifdef WITH_GUI
#include "mosaic_ui.hpp"
#endif

void usage() {
	std::cerr << "Usage: mosaic [options] -o <output file> <file [file [file [...]]]>\n"
		<< "    -o <file>        The output file.\n"
		<< "    -d <distance>    The feather distance in map units (default 100.)\n"
		<< "    <file [...]>     A list of files. The first is used as the background\n"
		<< "                     and determines the size of the output. Subsequent\n"
		<< "                     images are layered on top and feathered.\n"
		<< "    -v               Verbose messages.\n"
		<< "    -t               The number of threads to use. Defaults to the number\n"
		<< "                     of cores.\n"
		<< "    -s               Tile size. This manages the sizes of the tiles. Default 1024.\n"
		<< "                     that is used for feathering/blending. Keep in mind the width\n"
		<< "                     of the images, the number of threads and the size of the type.\n";
	
}

int runWithUI(int argc, char **argv) {
#ifdef WITH_GUI
	QApplication q(argc, argv);
	QWidget *w = new QWidget();
	geotools::ui::MosaicForm f;
	f.setupUi(w);
	w->show();
	return q.exec();
#else
	std::cerr << "GUI not enabled." << std::endl;
	return 1;
#endif
}

int main(int argc, char **argv) {

 	try {

	 	float distance = 100.0;
	 	std::vector<std::string> files;
	 	std::string outfile;
	 	int threads = 0;
	 	int tileSize = 1024;
	 	bool gui = false;

	 	g_loglevel(0);

	 	for(int i = 1; i < argc; ++i) {
	 		std::string arg(argv[i]);
	 		if(arg == "-d") {
	 			distance = atof(argv[++i]);
	 		} else if(arg == "-gui") {
	 			gui = true;
			} else if(arg == "-o") {
	 			outfile = argv[++i];
			} else if(arg == "-v") {
				g_loglevel(G_LOG_DEBUG);
			} else if(arg == "-t") {
				threads = atoi(argv[++i]);
			} else if(arg == "-s") {
				tileSize = atoi(argv[++i]);
	 		} else {
	 			files.push_back(argv[i]);
	 		}
	 	}

	 	if(gui) {
	 		return runWithUI(argc, argv);
	 	} else {
	 		geotools::raster::Mosaic m;
 			m.mosaic(files, outfile, distance, tileSize, threads);
 		}

 	} catch(const std::exception &e) {
 		g_error(e.what());
 		usage();
 		return 1;
 	}

 	return 0;
 }
