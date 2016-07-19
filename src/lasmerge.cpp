/*
 *
 *  Created on: Mar 13, 2015
 *      Author: rob
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <list>
#include <climits>
#include <memory>

#include "boost/filesystem.hpp"
#include "boost/algorithm/string/predicate.hpp"
#include "boost/algorithm/string/case_conv.hpp"

#include "liblas/liblas.hpp"

#include "geotools.h"

namespace las = liblas;

void usage() {
	std::cerr << "Usage: lasmerge [options] srcfiles*\n"
		<< " -o                              Speficy an output file.\n"
		<< " -v                              Verbose output.\n"
		<< " --maxx, --maxy, --minx, --miny  Specify the bounds of the data set.\n"
		<< "                                 Any combination of these can be used.\n";
}

void merge(std::vector<std::string> &files, std::string &outfile, double minx, double miny, double maxx, double maxy) {

	_trace("Processing " << files.size() << " files.");

	if(outfile.empty())
		_argerr("An output file (-o) is required.");

	if(files.size() == 0)
		_argerr("At least one input file is required.");

	if(minx >= maxx || miny >= maxy)
		_argerr("The bounds are silly.");

	las::Header *dsth = nullptr;
	las::Header::RecordsByReturnArray recs;
	las::ReaderFactory rf;

	// min x, max x, min y, max y, min z, max z
	double bounds[] = { DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_POS, DBL_MAX_NEG };
	int count = 0;
	std::vector<unsigned int> indices;

	for(unsigned int i = 0; i < files.size(); ++i) {

		_trace("Checking file " << files[i]);

		std::ifstream in(files[i], std::ios::in | std::ios::binary);
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();

		if(h.GetMinX() > maxx || h.GetMinY() > maxy || h.GetMaxX() < minx || h.GetMaxY() < miny)
			continue;
		
		if(h.GetMinX() < bounds[0]) bounds[0] = h.GetMinX();
		if(h.GetMaxX() > bounds[1]) bounds[1] = h.GetMaxX();
		if(h.GetMinY() < bounds[2]) bounds[2] = h.GetMinY();
		if(h.GetMaxY() > bounds[3]) bounds[3] = h.GetMaxY();
		if(h.GetMinZ() < bounds[4]) bounds[4] = h.GetMinZ();
		if(h.GetMaxZ() > bounds[5]) bounds[5] = h.GetMaxZ();

		indices.push_back(i);

		las::Header::RecordsByReturnArray rr = h.GetPointRecordsByReturnCount();

		if(dsth == nullptr) {
			// Create the destination header from the first source header. 
			// Initialize the return counts to zero.
			dsth = new las::Header(h);
			for(unsigned int i=0;i<rr.size();++i)
				recs.push_back(0);
		}

		// Rebuild the counts list because using the headers doesn't
		// always work.
		while(r.ReadNextPoint()) {
			las::Point pt = r.GetPoint();
			if(pt.GetX() > maxx || pt.GetX() < minx || pt.GetY() > maxy || pt.GetY() < miny)
				continue;
			++recs[pt.GetReturnNumber()-1];
			++count;
		}

		in.close();
	}

	if(indices.size() == 0)
		_argerr("No files matched the given bounds.")

	// Set the total count and update the point record counts.
	dsth->SetPointRecordsCount(count);
	for(unsigned int i=0;i<recs.size();++i)
		dsth->SetPointRecordsByReturnCount(i, recs[i]);

	dsth->SetMin(bounds[0], bounds[2], bounds[4]);
	dsth->SetMax(bounds[1], bounds[3], bounds[5]);

	std::ofstream out(outfile, std::ios::out | std::ios::binary);
	las::WriterFactory wf;
	las::Writer w(out, *dsth);

	_trace("Using points from " << indices.size() << " files.");

	for(unsigned int i = 0; i < indices.size(); ++i) {

		_trace("Processing file " << files[indices[i]]);

		std::string filename = files[indices[i]];
		std::ifstream in(filename, std::ios::in | std::ios::binary);
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();

		while(r.ReadNextPoint()) {
			las::Point pt = r.GetPoint();
			if(pt.GetX() > maxx || pt.GetX() < minx || pt.GetY() > maxy || pt.GetY() < miny)
				continue;
			w.WritePoint(r.GetPoint());
		}

		in.close();
	}

	out.close();
	delete dsth;	

}
int main(int argc, char ** argv) {

	try {
		std::vector<std::string> files;
		double minx = DBL_MAX_NEG;
		double maxx = DBL_MAX_POS;
		double miny = DBL_MAX_NEG;
		double maxy = DBL_MAX_POS;
		std::string outfile;
		bool verbose = false;

		for(int i=1;i<argc;++i) {
			std::string arg(argv[i]);
			if(arg == "--maxx") {
				maxx = atof(argv[++i]);
			} else if(arg == "--minx") {
				minx = atof(argv[++i]);
			} else if(arg == "--maxy") {
				maxy = atof(argv[++i]);
			} else if(arg == "--miny") {
				miny = atof(argv[++i]);
			} else if(arg == "-v") {
				verbose = true;
			} else if(arg == "-o") {
				outfile = argv[++i];
			} else {
				files.push_back(argv[i]);
			}
		}

		_loglevel(verbose ? LOG_TRACE : LOG_ERROR)

		merge(files, outfile, minx, miny, maxx, maxy);
	} catch(const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;

}
