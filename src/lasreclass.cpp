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
#include <string>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <liblas/liblas.hpp>

#include "geotools.h"

#define TIME_GAP 50.0

namespace las = liblas;

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

static int __seg_id = 0;
class Seg {
public:
	int id = ++__seg_id;
	double start = nan("");
	double end = nan("");
	las::Header *header = nullptr;
	las::Writer *writer = nullptr;
	std::ofstream out;
	bool inited = false;
	double bounds[6] = { DBL_MAX_POS, DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_NEG, DBL_MAX_POS, DBL_MAX_NEG };
	bool marked;
	int ptCount = 0;
	int ptRecByRetCount[5] = {0,0,0,0,0};
	Seg(double start, double end) {
		this->start = start;
		this->end = end;
	}

	void initialize(const std::string &filename, las::Header &h) {
		std::string file = filename + "_" + std::to_string(id) + ".las";
		if(!quiet)
			std::cerr << "Seg initialized as " << file << std::endl;
		out.open(file.c_str(), std::ios::out | std::ios::binary);
		header = new las::Header(h);
		writer = new las::Writer(out, *header);
		inited = true;
	}
	~Seg() {
		finalize();
	}
	void addPoint(las::Point &pt) {
		if(pt.GetX() < bounds[0]) bounds[0] = pt.GetX();
		if(pt.GetY() < bounds[1]) bounds[1] = pt.GetY();
		if(pt.GetX() > bounds[2]) bounds[2] = pt.GetX();
		if(pt.GetY() > bounds[3]) bounds[3] = pt.GetY();
		if(pt.GetZ() < bounds[4]) bounds[4] = pt.GetZ();
		if(pt.GetZ() > bounds[5]) bounds[5] = pt.GetZ();
		pt.SetPointSourceID(id);
		writer->WritePoint(pt);
		++ptCount;
		++ptRecByRetCount[pt.GetReturnNumber() - 1];
	}
	void finalize() {
		if(header != nullptr) {
			header->SetMin(bounds[0], bounds[1], bounds[4]);
			header->SetMax(bounds[2], bounds[3], bounds[5]);
			header->SetPointRecordsCount(ptCount);
			for(size_t i = 0; i < 5; ++i)
				header->SetPointRecordsByReturnCount(i + 1, ptRecByRetCount[i]);
			writer->SetHeader(*header);
			writer->WriteHeader();
			out.close();
			delete writer;
			delete header;
			header = nullptr;
			writer = nullptr;
		}
	}
	bool intersects(Seg *seg) {
		return !(seg->end < start || seg->start > end);
	}
	bool near(Seg *seg) {
		return ((start - seg->end) < 1.0 && start > seg->end) || ((seg->start - end) < 1.0 && seg->start > end);
	}
	bool insert(Seg *seg) {
		if(!quiet)
			std::cerr << "Insert: " << seg->id << ": " << seg->start << " > " << seg->end << std::endl;
		if(seg->id == id) {
			if(!quiet)
				std::cerr << "> same seg." << std::endl;
			return false;
		} else if(intersects(seg) || near(seg)) {
			if(!quiet)
				std::cerr << "> " << id << " joins " << seg->id << std::endl;
			start = _min(start, seg->start);
			end = _max(end, seg->end);
		} else {
			if(!quiet)
				std::cerr << "> " << id << " no relationship " << seg->id << std::endl;
			return false;
		}
		if(!quiet)
			std::cerr << "New span: " << id << ": " << start << " > " << end << std::endl;
		return true;
	}
};

struct {
  bool operator() (Seg *a, Seg *b) { 
  	return a->start < b->start;
  }
} segsort;

void normalizeFlightLines(std::vector<Seg*> &flightLines) {
	int count = flightLines.size();
	std::sort(flightLines.begin(), flightLines.end(), segsort);
	std::set<Seg*> output;
	Seg *seg = flightLines[0];
	output.insert(seg);
	for(size_t i = 1; i < flightLines.size(); ++i) {
		if(!seg->insert(flightLines[i])) {
			seg = flightLines[i];
			output.insert(seg);
		} else {
			delete flightLines[i];
		}
	}
	flightLines.resize(0);
	flightLines.assign(output.begin(), output.end());
	std::cerr << flightLines.size() << " remaining of " << count << std::endl;
	std::sort(flightLines.begin(), flightLines.end(), segsort);
	for(Seg *seg:flightLines)
		std::cerr << " -- " << seg->id << ": " << seg->start << " > " << seg->end << std::endl;
}

// TODO: Replace with interval tree.
int findFlightLine(std::vector<Seg*> &flightLines, double time) {
	for(Seg *seg:flightLines) {
		if(time >= seg->start && time <= seg->end)
			return seg->id;
	}
	if(!quiet)
		std::cerr << "Seg for time not found: " << time << std::endl;
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
			double time = pt.GetTime();
			if(startTime == 0) {
				startTime = endTime = time;
			} else {
				double gap = time - endTime;
				if(gap < 0.0 || gap > timeGap) { // One what?
					flightLines.push_back(new Seg(startTime, endTime));
					if(!quiet)
						std::cerr << "Time gap. Flightline: " << gap << "; " << startTime << " > " << endTime << std::endl;
					startTime = time;
				}
				endTime = time;
			}
  		}

		if(endTime != startTime) { // TODO: What is the time scale?)
			flightLines.push_back(new Seg(startTime, endTime));
			if(!quiet)
				std::cerr << "Time gap. Flightline: " << startTime << " > " << endTime << std::endl;
		}
  		in.close();
  	}

  	normalizeFlightLines(flightLines);


  	bool initSegs = true;
  	std::map<int, Seg*> flMap;

	for(unsigned int i = 0; i < files.size(); ++i) {

		std::string path = files[i];
		
		if(!quiet)
			std::cerr << "Processing " << path << std::endl;
		
		std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();

		if(initSegs) {
			std::string base(outfile + "/flight_line");
		  	for(Seg *seg:flightLines) {
		  		flMap[seg->id] = seg;
				seg->initialize(base, h);
		  	}
		  	initSegs = false;
		}
		
		while(r.ReadNextPoint()) {
			las::Point pt = r.GetPoint();
			Seg *seg = flMap[findFlightLine(flightLines, pt.GetTime())];
			seg->addPoint(pt);
  		}

		in.close();

	}

	for(Seg *seg:flightLines) {
		seg->finalize();
		delete seg;
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

		std::cerr << std::setprecision(12);

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
