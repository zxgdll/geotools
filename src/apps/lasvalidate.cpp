#include "lasvalidate.hpp"

#include "util.hpp"

void usage() {
	std::cerr << "Usage: lasvalidate [options] <lasfiles>\n"
			<< "    This program takes a spreadsheed of 3D survey locations and a LiDAR point cloud, in the form \n"
			<< "    of a list of LAS files, and interpolates the elevation of each survey point from the resulting \n"
			<< "    TIN. This allows the survey elevation to be compared to the LiDAR surface elevaton.\n"
			<< "    <lasfiles>  Is a list of las files.\n"
			<< "    -s  --  A csv file with x, y and z columns. These are the survey locations.\n"
			<< "    -o  --  The output file for survey points.\n"
			<< "    -p  --  The output file for LiDAR points (optional).\n"
			<< "    -r  --  The radius to search for lidar returns.\n"
			<< "    -c  --  The classes to include; comma-separated list.\n"
			<< "    -v  --  Verbose mode.\n";
}

int main(int argc, char **argv) {

 	try {

 		using namespace geotools::util;

 		std::vector<std::string> lasfiles;
  		std::string outfile;
  		std::string pointfile;
  		std::string datafile;
  		std::set<int> classes;
 		double distance = 0.0;

 		for(int i = 1; i < argc; ++i) {
 			std::string p(argv[i]);
 			if(p == "-o") {
 				outfile = argv[++i];
 			} else if(p == "-s") {
 				datafile = argv[++i];
 			} else if(p == "-p") {
 				pointfile = argv[++i];
 			} else if(p == "-r") {
 				distance = atof(argv[++i]);
 			} else if(p == "-v") {
 				g_loglevel(1);
 			} else if(p == "-c") {
				Util::intSplit(classes, argv[++i]);
 			} else {
 				lasfiles.push_back(argv[i]);
 			}
 		}

		geotools::las::validate(outfile, pointfile, datafile, lasfiles, classes, distance);

 	} catch(std::exception &e) {
 		std::cerr << e.what() << std::endl;
 		usage();
 		return 1;
 	}

 	return 0;
 }

