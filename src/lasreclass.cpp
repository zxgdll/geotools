/*
 * This tool reassigns the classes on each point.
 *
 *  Author: Rob Skelly rob@dijital.ca
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <list>
#include <climits>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <ogr_api.h>
#include <ogrsf_frmts.h>
#include <gdal_priv.h>

#include <geos/geom/GeometryCollection.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/Point.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CoordinateSequenceFactory.h>

#include <liblas/liblas.hpp>

namespace las = liblas;
namespace gg = geos::geom;

void usage() {
	std::cerr << "Usage: lasreclass [options] srcfiles*" << std::endl;
	std::cerr << "  Reassigns classes in each given las file." << std::endl;
	std::cerr << " -o Speficy an output folder (must exist)." << std::endl;
	std::cerr << " -m <from> <to> Specify a mapping. Points from the first class are reclassified to the second class. Repeatable." << std::endl;	
}

int main(int argc, char ** argv) {

	std::vector<std::string> files;
	std::map<int, int> mappings;
	std::string outfile;
	bool quiet = false;

	/* Parse and check input. */

	for(int i=1;i<argc;++i) {
		std::string arg(argv[i]);
		if(arg == "-m") {
			int from = atoi(argv[++i]);
			int to = atoi(argv[++i]);
			mappings[from] = to;
		} else if(arg == "-o") {
			outfile += argv[++i];
		} else if(arg == "-q") {
			quiet = true;
		} else {
			files.push_back(std::string(argv[i]));
		}
	}

	if(mappings.size() == 0) {
		std::cerr << "At least one mapping is required." << std::endl;
		usage();
		return 1;
	} else if(!quiet) {
		std::cerr << "Mappings: " << std::endl;
		for(std::map<int, int>::iterator it = mappings.begin(); it != mappings.end(); ++it)
			std::cerr << " " << it->first << " > " << it->second << std::endl;
	}

	if(outfile.size() == 0) {
		std::cerr << "An output file (-o) is required." << std::endl;
		usage();
		return 1;
	}

	if(files.size() == 0) {
		std::cerr << "At least one input file is required." << std::endl;
		usage();
		return 1;
	}

	/* Loop over files and figure out which ones are relevant. */
	las::ReaderFactory rf;

	for(unsigned int i = 0; i < files.size(); ++i) {

		std::string path = files[i];
		std::string base = path.substr(path.find_last_of("/") + 1);
		std::string newpath = outfile + "/" + base;
		
		if(!quiet)
			std::cerr << "Saving " << path << " to " << newpath << std::endl;
		
		std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();

		std::ofstream out(newpath.c_str(), std::ios::out | std::ios::binary);
		las::Header wh(h);
		las::Writer w(out, wh);


		while(r.ReadNextPoint()) {
			las::Point pt = r.GetPoint();
			las::Classification cl = pt.GetClassification();
			if(mappings.find(cl.GetClass()) != mappings.end()) {
				cl.SetClass(mappings[cl.GetClass()]);
				pt.SetClassification(cl);
			}
			w.WritePoint(pt);
		}

		in.close();
		out.close();
	}

	return 0;

}
