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
#include "RunningStats.hpp"

#define TIME_GAP 50.0

namespace las = liblas;
namespace gg = geos::geom;

void usage() {
	std::cerr << "Usage: lasreclass [options] srcfiles*\n"
				<< " Performs certain operations on LAS attributes.\n"
				<< " -o              Speficy an output folder (must exist).\n"
				<< " -m <from> <to>  Specify a class mapping. Points from the first class are reclassified to the second class. Repeatable.\n"
				<< " -f              Attempt to reconstruct flight lines using point times.\n";
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

void findFlightLine(std::map<double, int> &flightLines, double time, double *diff, int *id) {
	double nearTime = 0.0;
	int nearId = 0;
	for(auto it = flightLines.begin(); it != flightLines.end(); ++it) {
		if(std::abs(time - it->first) < nearTime) {
			nearTime = std::abs(time - it->first);
			nearId = it->second;
		}
	}
	*id = nearId;
	*diff = nearTime;
}

class Seg {
public:
	int id;
	double start = nan("");
	double end = nan("");
	Seg(int id, double start, double end) {
		this->id = id;
		this->start = start;
		this->end = end;
	}
	bool adjoins(Seg *seg) {
		return std::abs(end - seg->start) < 1.0 || std::abs(start - seg->end) < 1.0;
	}
	bool contains(Seg *seg) {
		return seg->start > start && seg->end < end;
	}
	bool insert(Seg *seg) {
		if(contains(seg)) {
			std::cerr << "> " << id << " contains " << seg->id << std::endl;
			return true;
		} else if(adjoins(seg)) {
			std::cerr << "> " << id << " adjoins " << seg->id << std::endl;
			if(seg->start >= end)
				end = seg->end;
			if(seg->end <= start)
				start = seg->start;
			return true;
		}
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
	for(int i = 0; i < flightLines.size(); ++i) {
		Seg *seg = flightLines[i];
		std::cerr << seg->id << ": " << seg->start << " > " << seg->end << std::endl;
		for(int j = i + 1; j < flightLines.size(); ++j) {
			if(seg->insert(flightLines[j])) {
				output.insert(seg);
			} else {
				i = j - 1;
				break;
			}
		}
	}
	flightLines.resize(0);
	flightLines.assign(output.begin(), output.end());
	std::cerr << flightLines.size() << " to remaining." << std::endl;
}

int findFlightLine(std::vector<Seg*> &flightLines, double time) {
	for(Seg *seg:flightLines) {
		if(time >= seg->start && time< seg->end)
			return seg->id;
	}
	return 0;
}

void recoverFlightlines(std::vector<std::string> &files, std::string &outfile, double stdDev) {

	if(outfile.size() == 0) 
		throw "An output directory (-o) is required.";

	if(files.size() == 0)
		throw "At least one input file is required.";

	if(stdDev <= 0)
		throw "Std. deviation must be larger than zero.";

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
				if(gap < 0.0 || gap > TIME_GAP) { // One what?
					flightLines.push_back(new Seg(pointSourceId, startTime, endTime));
					startTime = endTime;
					if(!quiet)
						std::cerr << "Time gap. Flightline: " << pointSourceId << ", " << gap << std::endl;
					++pointSourceId;
				}
			}
  		}

		if(endTime != startTime) { // One what?
			flightLines.push_back(new Seg(pointSourceId, startTime, endTime));
			if(!quiet)
				std::cerr << "Time gap. Flightline: " << pointSourceId << std::endl;
			++pointSourceId;
		}
  		in.close();
  	}

  	normalizeFlightLines(flightLines);

  	pointSourceId = 0;

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
			double time = pt.GetTime();
			pt.SetPointSourceID(findFlightLine(flightLines, time));
			w.WritePoint(pt);
  		}

		in.close();
		out.close();

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
	RunningStats rs;

	for(unsigned int i = 0; i < files.size(); ++i) {

		rs.clear();

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

	double flStdDev = 0;
	bool flightlines = false;
	
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
			flStdDev = atof(argv[++i]);
			flightlines = true;
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
			recoverFlightlines(files, outfile, flStdDev);
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
