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
#include <iomanip>
#include <queue>

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

#include "geotools.h"

#define TIME_GAP 50.0

namespace las = liblas;
namespace gg = geos::geom;

void usage() {
	std::cerr << "Usage: lasreclass [options] srcfiles*\n"
				<< " Performs certain operations on LAS attributes.\n"
				<< " -o              Speficy an output folder (must exist).\n"
				<< " -m <from> <to>  Specify a class mapping. Points from the first class are reclassified to the second class. Repeatable.\n"
				<< " -f              Attempt to reconstruct flight lines using point times.\n"
				<< " -fe             Export flightlines as single LAS files. Will have original filename with ID appended.\n"
				<< " -fg             The time gap required to signify a break in flight lines. Defaul " << TIME_GAP << "\n";
}

static bool quiet = true;

void mapClasses(std::vector<std::string> &files, std::string &outfile, std::map<int, int> &mappings) {

	if(mappings.size() == 0)
		throw "At least one mapping is required.";

	if(!quiet) {
		std::cerr << "Mappings: " << std::endl;
		for(std::map<int, int>::iterator it = mappings.begin(); it != mappings.end(); ++it)
			std::cerr << " " << it->first << " > " << it->second << std::endl;
	}

	if(outfile.size() == 0) 
		throw "An output directory (-o) is required.";

	if(files.size() == 0)
		throw "At least one input file is required.";

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
}

class Seg {
public:
	int id;
	double start = nan("");
	double end = nan("");
	las::Header *header;
	las::Writer *writer;
	std::ofstream out;
	bool inited = false;
	double bounds[6] = { DBL_MAX_POS, DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_NEG, DBL_MAX_POS, DBL_MAX_NEG };
	Seg(int id, double start, double end) {
		this->id = id;
		this->start = start;
		this->end = end;
	}

	void initialize(const std::string &filename, las::Header &h) {
		out.open(filename.c_str(), std::ios::out | std::ios::binary);
		header = new las::Header(h);
		writer = new las::Writer(out, *header);
		inited = true;
	}
	~Seg() {
		finalize();
	}
	void addPoint(las::Point &pt) {
		if(pt.GetX() < bounds[0]) bounds[0] = pt.GetX();
		if(pt.GetY() < bounds[1]) bounds[1] = pt.GetX();
		if(pt.GetX() > bounds[2]) bounds[2] = pt.GetX();
		if(pt.GetY() > bounds[3]) bounds[3] = pt.GetX();
		if(pt.GetZ() < bounds[4]) bounds[4] = pt.GetX();
		if(pt.GetZ() < bounds[5]) bounds[5] = pt.GetX();
		pt.SetPointSourceID(id);
		writer->WritePoint(pt);
	}
	void finalize() {
		if(header) {
			header->SetMin(bounds[0], bounds[1], bounds[4]);
			header->SetMax(bounds[1], bounds[3], bounds[5]);
			writer->WriteHeader();
			out.close();
			delete writer;
			delete header;
		}
	}
	bool adjoins(Seg *seg) {
		return std::abs(end - seg->start) < 1.0 || std::abs(start - seg->end) < 1.0;
	}
	bool contains(Seg *seg) {
		return seg->start > start && seg->end < end;
	}
	bool insert(Seg *seg) {
		if(contains(seg)) {
			if(!quiet)
				std::cerr << "> " << id << " contains " << seg->id << std::endl;
			return true;
		} else if(adjoins(seg)) {
			if(!quiet)
				std::cerr << "> " << id << " adjoins " << seg->id << std::endl;
			if(seg->start >= end)
				end = seg->end;
			if(seg->end <= start)
				start = seg->start;
			return true;
		}
		if(!quiet)
			std::cerr << "> " << id << " no relationship " << seg->id << std::endl;
		return false;
	}
};

struct {
  bool operator() (Seg *a, Seg *b) { 
  	return a->end < b->start;
  }
} segsort;

void normalizeFlightLines(std::vector<Seg*> &flightLines) {
	std::sort(flightLines.begin(), flightLines.end(), segsort);
	std::set<Seg*> output;
	for(size_t i = 0; i < flightLines.size(); ++i) {
		Seg *seg = flightLines[i];
		std::cerr << seg->id << ": " << seg->start << " > " << seg->end << std::endl;
		for(size_t j = i + 1; j < flightLines.size(); ++j) {
			if(seg->insert(flightLines[j])) {
				output.insert(seg);
			} else {
				delete flightLines[j];
				i = j - 1;
				break;
			}
		}
	}
	flightLines.resize(0);
	flightLines.assign(output.begin(), output.end());
	std::cerr << flightLines.size() << " to remaining." << std::endl;
}

// TODO: Replace with interval tree.
int findFlightLine(std::vector<Seg*> &flightLines, double time) {
	for(Seg *seg:flightLines) {
		if(time >= seg->start && time< seg->end)
			return seg->id;
	}
	return 0;
}

void recoverFlightlines(std::vector<std::string> &files, std::string &outfile, double timeGap, bool exportLas) {

	if(outfile.size() == 0) 
		throw "An output directory (-o) is required.";

	if(files.size() == 0)
		throw "At least one input file is required.";

	if(timeGap <= 0)
		throw "Time gap must be larger than zero.";

	/* Loop over files and figure out which ones are relevant. */
	las::ReaderFactory rf;

	// Store the IDs with the start or end time for that ID.
	// Subsequent files will search the table to find out if they have a 
	// continuation of this line.
	std::vector<Seg*> flightLines;
	int pointSourceId = 1;
	double startTime = 0, endTime = 0;

	for(unsigned int i = 0; i < files.size(); ++i) {
		std::string path = files[i];

		if(!quiet)
			std::cerr << "Checking " << path << std::endl;
		
		std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();

		startTime = 0;

		while(r.ReadNextPoint()) {
			las::Point pt = r.GetPoint();
			endTime = pt.GetTime();
			if(startTime == 0) {
				startTime = endTime;
			} else {
				double gap = endTime - startTime;
				if(gap < 0.0 || gap > timeGap) { // One what?
					flightLines.push_back(new Seg(pointSourceId, startTime, endTime));
					startTime = endTime;
					if(!quiet)
						std::cerr << "Time gap. Flightline: " << pointSourceId << ", " << gap << std::endl;
					++pointSourceId;
				}
			}
  		}

		if(endTime != startTime) { // TODO: What is the time scale?
			flightLines.push_back(new Seg(pointSourceId, startTime, endTime));
			if(!quiet)
				std::cerr << "Time gap. Flightline: " << pointSourceId << std::endl;
			++pointSourceId;
		}
  		in.close();
  	}

  	normalizeFlightLines(flightLines);

  	std::map<int, Seg*> flMap;
  	for(Seg *seg:flightLines)
  		flMap[seg->id] = seg;

  	pointSourceId = 0;

	for(unsigned int i = 0; i < files.size(); ++i) {

		std::string path = files[i];
		std::string base(outfile + "/" + path.substr(path.find_last_of("/") + 1, path.find_last_of(".")));
		
		if(!quiet)
			std::cerr << "Saving " << path << " to " << base << std::endl;
		
		std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();

		//std::ofstream out((base + ".las").c_str(), std::ios::out | std::ios::binary);
		//las::Header wh(h);
		//las::Writer w(out, wh);
		
		while(r.ReadNextPoint()) {
			las::Point pt = r.GetPoint();
			double time = pt.GetTime();
			Seg *seg = flMap[findFlightLine(flightLines, time)];
			if(!seg->inited)
				seg->initialize(base, h);
			seg->addPoint(pt);
  		}

  		for(Seg *seg:flightLines) {
  			seg->finalize();
  			delete seg;
  		}

		in.close();

	}
}

void computeDirection(std::queue<las::Point> &pq, double *qdir, double *qvel) {
	bool first = true;
	double dir = 0;
	double vel = 0;
	las::Point *pt0;
	std::queue<las::Point> pq0(pq);
	while(pq0.size()) {
		las::Point &pt = pq0.front();
		pq0.pop();
		if(first) {
			first = false;
		} else {
			dir += atan2(pt.GetY() - pt0->GetY(), pt.GetX() - pt0->GetX());
			vel += sqrt(_sq(pt.GetY() - pt0->GetY()) + _sq(pt.GetX() - pt0->GetX()));
		}
		pt0 = &pt;
	}
	*qdir = dir / 5.0;
	*qvel = vel / 5.0;
}

void computeDirection(las::Point &pt0, las::Point &pt, double *pdir, double *pvel) {
	*pdir = atan2(pt.GetY() - pt0.GetY(), pt.GetX() - pt0.GetX());
	*pvel = sqrt(_sq(pt.GetY() - pt0.GetY()) + _sq(pt.GetX() - pt0.GetX()));
}

void recoverEdges(std::vector<std::string> &files, std::string &outfile) {

	if(outfile.size() == 0) 
		throw "An output directory (-o) is required.";

	if(files.size() == 0)
		throw "At least one input file is required.";

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

		//int lastPSId = 0;
		double qdir = 0, qvel = 0; // queue direction
		double pdir = 0, pvel = 0; // point direction
		std::queue<las::Point> pq;

		while(r.ReadNextPoint()) {
			las::Point pt(r.GetPoint());
			if(pq.size() == 5) {
				computeDirection(pq.front(), pq.back(), &qdir, &qvel);
				computeDirection(pq.back(), pt, &pdir, &pvel);
				pq.pop();
				while(qdir < 0)
					qdir += PI * 2;
				while(pdir < 0)
					pdir += PI * 2;
				//std::cerr << _deg(qdir) << ", " << qvel << "; " << _deg(pdir) << ", " << pvel << std::endl;
				if(std::abs(qdir - pdir) > PI * 0.75)
					pt.GetClassification().SetClass(31);
			}
			pq.push(pt);
			w.WritePoint(pt);
  		}

		in.close();
		out.close();

	}
}

int main(int argc, char ** argv) {

	std::vector<std::string> files;
	std::map<int, int> mappings;
	bool mapping = false;

	bool flightlines = false;
	double flTimeGap = TIME_GAP;
	bool flExport = false;

	bool edges = false;

	std::string outfile;

	/* Parse and check input. */

	for(int i=1;i<argc;++i) {
		std::string arg(argv[i]);
		if(arg == "-m") {
			int from = atoi(argv[++i]);
			int to = atoi(argv[++i]);
			mappings[from] = to;
			mapping = true;
		} else if(arg == "-f") {
			flightlines = true;
		} else if(arg == "-fg") {
			flTimeGap = atof(argv[++i]);
		} else if(arg == "-fe") {
			flExport = true;
		} else if(arg == "-o") {
			outfile += argv[++i];
		} else if(arg == "-e") {
			edges = true;
		} else if(arg == "-v") {
			quiet = false;
		} else {
			files.push_back(std::string(argv[i]));
		}
	}

	try {
		if(mapping) {
			mapClasses(files, outfile, mappings);
		} else if(flightlines) {
			recoverFlightlines(files, outfile, flTimeGap, flExport);
		} else if(edges) {
			recoverEdges(files, outfile);
		} else {
			throw "No command given.";
		}
	} catch(const char *e) {
		std::cerr << e << std::endl;
		usage();
		return 1;
	}

	return 0;

}
