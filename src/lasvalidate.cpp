#include <liblas/liblas.hpp>

#include "geotools.h"
#include "Util.hpp"

namespace las = liblas;

void usage() {
	std::cerr << "Usage: lasvalidate [options] <lasfiles>" << std::endl
			<< "    <lasfiles>  Is a list of las files." << std::endl
			<< "    -p  --  A csv file with x, y and z columns. These are the survey locations." << std::endl
			<< "    -o  --  The output file." << std::endl
			<< "    -r  --  The radius to search for lidar returns (not compatible with -n)." << std::endl
			<< "    -n  --  The number of neighbours (not compatible with -r)." << std::endl
			<< "    -c  --  The classes to include; comma-separated list." << std::endl;
}

double dist(double x1, double y1, double x2, double y2) {
	return sqrt(_sq(x2-x1) + _sq(y2-y1));
}

void validate(std::string &outfile, std::string &datafile, std::vector<std::string> &lasfiles, 
	std::set<int> &classes, double distance, int neighbours) {

	if(outfile.empty())
		throw "Outfile not given.";
	if(datafile.empty())
		throw "Data file not given.";
	if(lasfiles.size() == 0)
		throw "No las files.";
	if(distance != 0.0 && neighbours != 0)
		throw "Distance and neighbours arguments are mutually exclusive.";
	if(distance < 0.0)
		throw "Distance must be greater than zero.";
	if(neighbours < 0)
		throw "Neighbours must be greater than zero.";

	std::vector<std::tuple<double, double, double> > points;
	Util::loadXYZSamples(datafile, points);

	las::ReaderFactory rf;

	for(unsigned int i=0; i<lasfiles.size(); ++i) {
		std::cout << lasfiles[i] << std::endl;
		std::ifstream in(lasfiles[i].c_str());
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();

		while(r.ReadNextPoint()) {
			las::Point pt = r.GetPoint();
			// If this point is not in the class list, skip it.
			int cls = pt.GetClassification().GetClass();
			if(classes.size() > 0 && !Util::inList(classes, cls))
				continue;

			double lx = pt.GetX();
			double ly = pt.GetY();

			for(unsigned int j = 0; j < points.size(); ++j) {
				double px = std::get<0>(points[i]);
				double py = std::get<1>(points[i]);
				// If outside of the radius, skip it.
				if(dist(px, py, lx, ly) > distance)
					continue;

				double pz = std::get<2>(points[i]);
				double lz = pt.GetZ();

				// Output
				std::cout << j << "," << px << "," << py << "," << pz << "," << lx << "," << ly << "," << lz << "," << cls << std::endl;

			}
		}
	}


}

int main(int argc, char **argv) {

 	try {

 		std::vector<std::string> lasfiles;
  		std::string outfile;
  		std::string datafile;
  		std::set<int> classes;
 		double distance = 0.0;
 		int neighbours = 0;

 		for(int i = 0; i < argc; ++i) {
 			std::string p(argv[i]);
 			if(p == "-o") {
 				outfile = argv[++i];
 			} else if(p == "-p") {
 				datafile = argv[++i];
 			} else if(p == "-r") {
 				distance = atof(argv[++i]);
 			} else if(p == "-n") {
 				neighbours = atoi(argv[++i]);
 			} else if(p == "-c") {
				Util::intSplit(classes, argv[++i]);
 			} else {
 				lasfiles.push_back(argv[i]);
 			}
 		}

		validate(outfile, datafile, lasfiles, classes, distance, neighbours);

 	} catch(const char *e) {
 		std::cerr << e << std::endl;
 		usage();
 		return 1;
 	}

 	return 0;
 }



