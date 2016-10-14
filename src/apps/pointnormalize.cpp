#include <iostream>

#include "geotools.h"
#include "pointnormalize.hpp"

using namespace geotools::point;

void usage() {
	std::cerr << "Usage: pointnormalize <options> <terrain file> <point file [point file [point file ...]]>\n"
		<< " --point-dir                 The directory where normalized point files are to be written.\n"
		<< " -v                          Verbose output.\n"
		<< " -h                          Print this message.\n"
		<< " --threads                   The number of threads to use for computing output.\n"
		<< " -gui                        Run the graphical user interface.\n";
}
	
int main(int argc, char **argv) {

	try {

		std::string terrainFile;
		std::string chmFile;
		std::string pointDir;
		std::list<std::string> pointFiles;
		int threads = 1;
		bool gui = false;

		g_loglevel(0);
		
		for(int i = 1; i < argc; ++i) {
			std::string s(argv[i]);
			if(s == "-h") {
				usage();
				return 0;
			} else if(s == "-gui") {
				gui = true;
			} else if(s == "-v") {
				g_loglevel(G_LOG_DEBUG);
			} else if(s == "--threads") {
				threads = atoi(argv[++i]);
			} else if(s == "--point-dir") {
				pointDir.assign(argv[++i]);
			} else {
				if(terrainFile.empty()) {
					terrainFile.assign(s);
				} else {
					pointFiles.push_back(s);
				}
			}
		}

		if(gui) {
			//return runWithUI(argc, argv);
		} else {
			PointNormalize pn;
			PointNormalizeConfig config;
			config.terrainFile = terrainFile;
			config.pointOutputDir = pointDir;
			config.pointFiles = pointFiles;
			config.threads = threads;
			
			pn.normalize(config);
		}

	} catch(const std::exception &ex) {
		g_error(ex.what());
		usage();
		return 1;
	}

	return 0;
}
