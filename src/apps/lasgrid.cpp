
#include <set>
#include <vector>

#include "lasgrid.hpp"
#include "lasgrid_ui.hpp"

using namespace geotools::util;
using namespace geotools::las;
using namespace geotools::las::lasgrid_util;

void usage() {
	std::cerr << "Usage: lasgrid <options> <file [file [file]]>\n"
		<< " -o <output file>\n"
		<< " -t <type>                   Output median, mean, max, min, variance (sample), pvariance (population),\n"
		<< "                             count, density, stddev (sample), pstddev (population). Default mean.\n"
		<< " -r <resolution>             Resolution (default 2).\n"
		<< " -s <srid>                   The EPSG ID of the CRS.\n"
		<< " -c <classes>                Comma-delimited (e.g. '2,0' (ground and unclassified)).\n"
		<< " -a <attribute>              Use height, intensity (default height).\n"
		<< " -d <radius>                 Radius (not diameter); use zero for cell bounds.\n"
		<< "                             For example, if the cell size is 2, the circumcircle's radius is sqrt(2) (~1.41).\n"
		<< " -b <minx miny maxx maxy>    Extract points from the given box and create a raster of this size.\n"
		<< " -f                          Fill voids.\n"
		<< " -v                          Verbose output.\n"
		<< " -h                          Print this message.\n"
		<< " --angle-limit               Points located outside of this angle (devation from nadir) are excluded.\n"
		<< " -gui                        Run the graphical user interface.\n";
}

int runWithUI(int argc, char **argv) {
	QApplication q(argc, argv);
	QWidget *w = new QWidget();
	geotools::ui::LasgridForm f;
	f.setupUi(w);
	w->show();
	return q.exec();
}

int main(int argc, char **argv) {

	try {

		std::string dstFile;
		int crs = 0;
		int type = TYPE_MEAN;
		int att = ATT_HEIGHT;
		bool fill = false;
		double resolution = 2.0;
		double radius = -1.0;
		unsigned char angleLimit = 100;
		Bounds bounds;
		std::set<int> classes;
		std::vector<std::string> files;

		for(int i = 1; i < argc; ++i) {
			std::string s(argv[i]);
			if(s == "-h") {
				usage();
				return 0;
			} else if(s == "-gui") {
				return runWithUI(argc, argv);
			} else if(s == "-o") {
				dstFile = argv[++i];
			} else if(s == "-s") {
				crs = atoi(argv[++i]);
			} else if(s == "-f") {
				fill = true;
			} else if(s == "-t") {
				type = parseType(argv[++i]);
			} else if(s == "-r") {
				resolution = atof(argv[++i]);
			} else if(s == "-c") {
				Util::intSplit(classes, argv[++i]);
			} else if(s == "-a") {
				att = parseAtt(argv[++i]);
			} else if(s == "-d") {
				radius = atof(argv[++i]);
			} else if(s == "-v") {
				g_loglevel(G_LOG_TRACE);
			} else if(s == "--angle-limit") {
				angleLimit = (unsigned char) atoi(argv[++i]);
			} else if(s == "-b") {
				bounds.extend(atof(argv[++i]), atof(argv[++i]));
				bounds.extend(atof(argv[++i]), atof(argv[++i]));
			} else {
				files.push_back(argv[i]);
			}
		}

		lasgrid(dstFile, files, classes, crs, att, type, radius, resolution, bounds, angleLimit, fill);

	} catch(const std::exception &ex) {
		std::cerr << ex.what() << std::endl;
		usage();
		return 1;
	}

	return 0;
}
