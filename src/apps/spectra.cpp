#include <iostream>

#include "omp.h"

#include "spectra.hpp"

using namespace geotools::spectral;
using namespace geotools::spectral::config;

void usage() {
	std::cerr << "Usage" << std::endl;
}

int main(int argc, char **argv) {

	try {

		SpectralConfig config;		
		int threads = 1;

		for(int i = 1; i < argc; ++i) {
			std::string arg(argv[i]);
			if(arg == "-i") {
				config.indexFilename.assign(argv[++i]);
			} else if(arg == "-o") {
				config.outputFilename.assign(argv[++i]);
			} else if(arg == "-b") {
				config.bands.insert(atoi(argv[++i]));
			} else if(arg == "-n") {
				config.nodata = atof(argv[++i]);
			} else if(arg == "-threads") {
				threads = atoi(argv[++i]);
			} else {
				config.spectralFilenames.push_back(std::string(argv[i]));
			}	
		}

		if(threads <= 0) 
			threads = 1;

		omp_set_dynamic(0);
		omp_set_num_threads(threads);

		Spectral s;
		s.extractSpectra(config);

	} catch(const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;
}
