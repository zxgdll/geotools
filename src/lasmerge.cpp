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


namespace las = liblas;

void usage() {
	std::cerr << "Usage: lasmerge [options] srcfiles*" << std::endl;
	std::cerr << " -o Speficy an output file." << std::endl;
	std::cerr << " --maxx, --maxy, -- minx, --miny Specify the bounds of the data set. Any combination of these can be used." << std::endl;	
}

int main(int argc, char ** argv) {

	std::vector<char *> files;
	double minx = DBL_MAX_NEG;
	double maxx = DBL_MAX_POS;
	double miny = DBL_MAX_NEG;
	double maxy = DBL_MAX_POS;
	char *outfile = NULL;
	bool quiet = false;

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
		} else if(arg == "-q") {
			quiet = true;
		} else if(arg == "-o") {
			outfile = argv[++i];
		} else {
			files.push_back(argv[i]);
		}
	}

	if(!quiet)
		std::cerr << "Processing " << files.size() << " files." << std::endl;

	if(outfile == NULL) {
		std::cerr << "An output file (-o) is required." << std::endl;
		usage();
		return 1;
	}

	if(files.size() == 0) {
		std::cerr << "At least one input file is required." << std::endl;
		usage();
		return 1;
	}

	if(minx >= maxx || miny >= maxy) {
		std::cerr << "The bounds are silly: " << minx << ", " << miny << ", " << maxx << ", " << maxy << std::endl;
		usage();
		return 1;
	}

	las::Header * dsth = NULL;
	las::Header::RecordsByReturnArray recs;
	las::ReaderFactory rf;

	// min x, max x, min y, max y, min z, max z
	double bounds[] = { DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_POS, DBL_MAX_NEG };
	int count = 0;
	std::vector<unsigned int> indices;

	for(unsigned int i = 0; i < files.size(); ++i) {

		const char * filename = files[i];
		std::ifstream in(filename, std::ios::in | std::ios::binary);
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();

		if(!quiet)
			std::cout << "Checking file " << filename << std::endl;

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

		if(dsth == NULL) {
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

	if(indices.size() == 0) {
		std::cerr << "No files matched the given bounds." << std::endl;
		return 1;
	}

	// Set the total count and update the point record counts.
	dsth->SetPointRecordsCount(count);
	for(unsigned int i=0;i<recs.size();++i)
		dsth->SetPointRecordsByReturnCount(i, recs[i]);

	dsth->SetMin(bounds[0], bounds[2], bounds[4]);
	dsth->SetMax(bounds[1], bounds[3], bounds[5]);

	std::ofstream out(outfile, std::ios::out | std::ios::binary);
	las::WriterFactory wf;
	las::Writer w(out, *dsth);

	if(!quiet)
		std::cout << "Using points from " << indices.size() << " files." << std::endl;
	for(unsigned int i = 0; i < indices.size(); ++i) {

		const char * filename = files[indices[i]];
		std::ifstream in(filename, std::ios::in | std::ios::binary);
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();

		if(!quiet)
			std::cout << "Processing file " << filename << std::endl;

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

	return 0;

}
