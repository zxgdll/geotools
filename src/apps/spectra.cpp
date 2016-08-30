#include <iostream>

#include "omp.h"

#include "spectra.hpp"

using namespace geotools::spectral;
using namespace geotools::spectral::config;

void usage() {
	std::cerr << "This program extracts digital numbers from a configurable range\n"
		<< "of bands in one raster, grouped by the IDs stored in another raster.\n"
		<< "the IDs in the ID raster are assumed to represent contiguous polygons.\n"
		<< "The intended usage is to extract spectra from tree crowns.\n\n"
		<< "Usage: spectra-app <options> <spectral file [spectral file [...]]>\n"
		<< " -b  -- Band configuration. A comma-separated list of individual bands,\n"
		<< "        or ranges of bands indicated by a dash.\n"
		<< "        example: 1,2,5-22,100 would be bands 1, 2 and 100 plus all bands\n"
		<< "        if not given, all bands are used.\n"
		<< "        between 5 and 22 (inclusive.)\n"
		<< " -i  -- ID raster. Unsigned 32 bit integer.\n"
		<< " -in -- ID raster nodata. Default 0.\n"
		<< " -o  -- Output file; a SpatiaLite database containing the pixel-center\n"
		<< "        geometry, the ID and a column for each band populated with the \n"
		<< "        digital number of the spectral band at that location.\n"
		<< " -sn -- The spectral nodata value. Expected to be appropriate for the\n"
		<< "        spectral raster's data type.\n";
}

int main(int argc, char **argv) {

	try {

		SpectralConfig config;		
		int threads = 0;

		for(int i = 1; i < argc; ++i) {
			std::string arg(argv[i]);
			if(arg == "-i") {
				config.indexFilename.assign(argv[++i]);
			} else if(arg == "-o") {
				config.outputFilename.assign(argv[++i]);
			} else if(arg == "-b") {
				config.bands.insert(atoi(argv[++i]));
			} else if(arg == "-sn") {
				config.nodata = atof(argv[++i]);
			} else if(arg == "-threads") {
				threads = atoi(argv[++i]);
			} else {
				config.spectralFilenames.push_back(std::string(argv[i]));
			}	
		}

		if(threads <= 0) 
			threads = 0;

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
