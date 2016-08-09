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

namespace geotools {

	namespace las {

		void merge(std::vector<std::string> &files, std::string &outfile, double minx, double miny, double maxx, double maxy) {

			g_trace("Processing " << files.size() << " files.");

			if(outfile.empty())
				g_argerr("An output file (-o) is required.");

			if(files.size() == 0)
				g_argerr("At least one input file is required.");

			if(minx >= maxx || miny >= maxy)
				g_argerr("The bounds are silly.");

			liblas::Header *dsth = nullptr;
			liblas::Header::RecordsByReturnArray recs;
			liblas::ReaderFactory rf;

			// min x, max x, min y, max y, min z, max z
			double bounds[] = { G_DBL_MAX_POS, G_DBL_MAX_NEG, G_DBL_MAX_POS, G_DBL_MAX_NEG, G_DBL_MAX_POS, G_DBL_MAX_NEG };
			int count = 0;
			std::vector<unsigned int> indices;

			for(unsigned int i = 0; i < files.size(); ++i) {

				g_trace("Checking file " << files[i]);

				std::ifstream in(files[i], std::ios::in | std::ios::binary);
				liblas::Reader r = rf.CreateWithStream(in);
				liblas::Header h = r.GetHeader();

				if(h.GetMinX() > maxx || h.GetMinY() > maxy || h.GetMaxX() < minx || h.GetMaxY() < miny)
					continue;
				
				if(h.GetMinX() < bounds[0]) bounds[0] = h.GetMinX();
				if(h.GetMaxX() > bounds[1]) bounds[1] = h.GetMaxX();
				if(h.GetMinY() < bounds[2]) bounds[2] = h.GetMinY();
				if(h.GetMaxY() > bounds[3]) bounds[3] = h.GetMaxY();
				if(h.GetMinZ() < bounds[4]) bounds[4] = h.GetMinZ();
				if(h.GetMaxZ() > bounds[5]) bounds[5] = h.GetMaxZ();

				indices.push_back(i);

				liblas::Header::RecordsByReturnArray rr = h.GetPointRecordsByReturnCount();

				if(dsth == nullptr) {
					// Create the destination header from the first source header. 
					// Initialize the return counts to zero.
					dsth = new liblas::Header(h);
					for(unsigned int i=0;i<rr.size();++i)
						recs.push_back(0);
				}

				// Rebuild the counts list because using the headers doesn't
				// always work.
				while(r.ReadNextPoint()) {
					liblas::Point pt = r.GetPoint();
					if(pt.GetX() > maxx || pt.GetX() < minx || pt.GetY() > maxy || pt.GetY() < miny)
						continue;
					++recs[pt.GetReturnNumber()-1];
					++count;
				}

				in.close();
			}

			if(indices.size() == 0)
				g_argerr("No files matched the given bounds.")

			// Set the total count and update the point record counts.
			dsth->SetPointRecordsCount(count);
			for(unsigned int i=0;i<recs.size();++i)
				dsth->SetPointRecordsByReturnCount(i, recs[i]);

			dsth->SetMin(bounds[0], bounds[2], bounds[4]);
			dsth->SetMax(bounds[1], bounds[3], bounds[5]);

			std::ofstream out(outfile, std::ios::out | std::ios::binary);
			liblas::WriterFactory wf;
			liblas::Writer w(out, *dsth);

			g_trace("Using points from " << indices.size() << " files.");

			for(unsigned int i = 0; i < indices.size(); ++i) {

				g_trace("Processing file " << files[indices[i]]);

				std::string filename = files[indices[i]];
				std::ifstream in(filename, std::ios::in | std::ios::binary);
				liblas::Reader r = rf.CreateWithStream(in);
				liblas::Header h = r.GetHeader();

				while(r.ReadNextPoint()) {
					liblas::Point pt = r.GetPoint();
					if(pt.GetX() > maxx || pt.GetX() < minx || pt.GetY() > maxy || pt.GetY() < miny)
						continue;
					w.WritePoint(r.GetPoint());
				}

				in.close();
			}

			out.close();
			delete dsth;	

		}

	} // las

} // geotools



void usage() {
	std::cerr << "Usage: lasmerge [options] srcfiles*\n"
		<< " -o                              Speficy an output file.\n"
		<< " -v                              Verbose output.\n"
		<< " --maxx, --maxy, --minx, --miny  Specify the bounds of the data set.\n"
		<< "                                 Any combination of these can be used.\n";
}


int main(int argc, char ** argv) {

	std::vector<std::string> files;
	double minx = G_DBL_MAX_NEG;
	double maxx = G_DBL_MAX_POS;
	double miny = G_DBL_MAX_NEG;
	double maxy = G_DBL_MAX_POS;
	std::string outfile;

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
			g_loglevel(1);
		} else if(arg == "-o") {
			outfile = argv[++i];
		} else {
			files.push_back(argv[i]);
		}
	}

	try {

		geotools::las::merge(files, outfile, minx, miny, maxx, maxy);
		
	} catch(const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;

}
