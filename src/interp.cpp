/*
 *  Created on: Feb 22, 2016
 *  Author: rob
 */

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <algorithm>

#include "csv.h"

#include "Util.hpp"
#include "Raster.hpp"
#include "interp/InterpPoint.hpp"
#include "interp/Interpolator.hpp"
#include "interp/IDWInterpolator.hpp"
#include "interp/AvgInterpolator.hpp"
#include "interp/PlanarInterpolator.hpp"
#include "interp/NaturalNeighbourInterpolator.hpp"
#include "interp/SimpleKrigingInterpolator.hpp"

//#include "ShapeWriter.hpp"

//ShapeWriter sw("/home/rob/Documents/geotools/data/out.shp");

double _min(double a, double b) {
	return a > b ? b : a;
}

double _max(double a, double b) {
	return a < b ? b : a;
}

double _random() {
	double r = ((double) std::rand()) / RAND_MAX;
	return r;
}

double _sq(double a) {
	return a*a;
}

void loadSamples(std::string &datafile, std::list<interp::InterpPoint> &samples) {
	io::CSVReader<3> in(datafile.c_str());
	in.read_header(io::ignore_extra_column, "x", "y", "z");
	double x, y, z;
	while(in.read_row(x, y, z))
		samples.push_back(interp::InterpPoint(x, y, z));
}

void interpolate(std::string &datafile, std::string &templatefile, std::string &outfile,
		interp::Interpolator &inter,
		double resolution, bool print) {

	if(templatefile == outfile)
		throw "The output file must not be the same as the template file.";

	if(resolution <= 0.0)
		throw "Invalid resolution.";

	// Initializes source, destination rasters.
	Raster<float> tpl(templatefile, 1, false);
	std::string proj;
	tpl.projection(proj);
	Raster<float> out(outfile, tpl.minx(), tpl.miny(), tpl.maxx(), tpl.maxy(),
			resolution, tpl.nodata(), proj.c_str());

	std::list<interp::InterpPoint> samples;

	loadSamples(datafile, samples);

	if(samples.size() == 0)
		throw "No samples loaded.";

	std::cerr << "Interpolating..." << std::endl;
	inter.interpolate(out, samples);

}
void usage() {
	std::cerr << "Usage: interp [options]" << std::endl
			<< " -t -- type:          The type of adjustment. " << std::endl
			<< "                      nn  - natural neighbours." << std::endl
			<< "                      pl  - plane fit. " << std::endl
			<< "                      avg - shift vertically by the average difference. " << std::endl
			<< "                      idw - inverse distance weighting (use -e switch for exponent; default 1). " << std::endl
			<< "                      sk  - Simple Kriging." << std::endl
			<< " -r -- resolution:    the pixel size of the output." << std::endl
			<< " -i -- template file: A template file to produce the output file." << std::endl
			<< " -o -- output file:   The output file." << std::endl
			<< " -d -- data file:     A CSV file with data." << std::endl
			<< " -e -- idw exponent." << std::endl;
}

int main(int argc, char **argv) {

 	try {

 		std::string templatefile;
  		std::string outfile;
  		std::string datafile;
  		std::string type;
 		double resolution = 0.0;
 		bool print = false;
 		double idwExp = 1.0;

 		for(int i = 0; i < argc; ++i) {
 			std::string p(argv[i]);
 			if(p == "-t") {
 				type = argv[++i];
 			} else if(p == "-o") {
 				outfile = argv[++i];
 			} else if(p == "-i") {
 				templatefile = argv[++i];
 			} else if(p == "-r") {
 				resolution = atof(argv[++i]);
 			} else if(p == "-p") {
 				print = true;
 			} else if(p == "-e") {
 				idwExp = atof(argv[++i]);
 			} else if(p == "-d") {
 				datafile = argv[++i];
 			}
 		}

 		if(resolution <= 0)
 			throw "Invalid resolution.";
 		if(outfile.empty())
 			throw "An output file is required.";
 		if(templatefile.empty())
 			throw "A template file is required.";
 		if(datafile.empty())
 			throw "A data file is required.";

 		interp::Interpolator *inter;

 		if(type == "nn") {
 			inter = new interp::naturalneighbour::NaturalNeighbourInterpolator();
 		} else if(type == "pl") {
 			inter = new interp::planar::PlanarInterpolator();
 		} else if(type == "avg") {
 			inter = new interp::avg::AvgInterpolator();
 		} else if(type == "idw") {
 			inter = new interp::idw::IDWInterpolator(idwExp);
 		} else if(type == "sk") {
 			interp::kriging::SimpleKrigingInterpolator *sk = new interp::kriging::SimpleKrigingInterpolator(argc, argv);
 			inter = sk;
 		} else {
 			throw "Unknown interpolator type.";
 		}

		interpolate(datafile, templatefile, outfile, *inter, resolution, print);

		delete inter;

 		std::cerr << "Done." << std::endl;

 	} catch(const char *e) {
 		std::cerr << e << std::endl;
 		usage();
 		return 1;
 	}

 	return 0;
 }



