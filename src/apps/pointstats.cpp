
#include <set>
#include <vector>

#include "pointstats.hpp"
#include "pointstats_ui.hpp"

using namespace geotools::util;
using namespace geotools::point;

void usage() {
	std::cerr << "Usage: pointstats <options> <file [file [file]]>\n"
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
	geotools::ui::PointStatsForm f;
	f.setupUi(w);
	w->show();
	return q.exec();
}

int main(int argc, char **argv) {

	try {

		std::string dstFile;
		int crs = 0;
		std::string type = "mean";
		std::string att = "height";
		bool fill = false;
		double resolution = 2.0;
		double radius = -1.0;
		unsigned char angleLimit = 100;
		Bounds bounds;
		std::set<unsigned char> classes;
		std::list<std::string> files;
		bool gui = false;

		g_loglevel(0);
		
		for(int i = 1; i < argc; ++i) {
			std::string s(argv[i]);
			if(s == "-h") {
				usage();
				return 0;
			} else if(s == "-gui") {
				gui = true;
			} else if(s == "-o") {
				dstFile = argv[++i];
			} else if(s == "-s") {
				crs = atoi(argv[++i]);
			} else if(s == "-f") {
				fill = true;
			} else if(s == "-t") {
				type = argv[++i];
			} else if(s == "-r") {
				resolution = atof(argv[++i]);
			} else if(s == "-c") {
				Util::intSplit(classes, argv[++i]);
			} else if(s == "-a") {
				att = argv[++i];
			} else if(s == "-d") {
				radius = atof(argv[++i]);
			} else if(s == "-v") {
				g_loglevel(G_LOG_DEBUG);
			} else if(s == "--angle-limit") {
				angleLimit = (unsigned char) atoi(argv[++i]);
			} else if(s == "-b") {
				bounds.extend(atof(argv[i + 1]), atof(argv[i + 2]));
				bounds.extend(atof(argv[i + 3]), atof(argv[i + 4]));
				i += 4;
			} else {
				files.push_back(argv[i]);
			}
		}

		if(gui) {
			return runWithUI(argc, argv);
		} else {
			PointStats lg;
			PointStatsConfig config;
			config.dstFile = dstFile;
			config.sourceFiles = files;
			config.classes = classes;
			config.hsrid = crs;
			config.attribute = config.parseAtt(att);
			config.type = config.parseType(type);
			config.radius = radius;
			config.resolution = resolution;
			config.bounds = bounds;
			config.angleLimit = angleLimit;
			config.fill = fill;
			config.snap = true; // TODO: Snap
			lg.pointstats(config);
		}

	} catch(const std::exception &ex) {
		g_error(ex.what());
		usage();
		return 1;
	}

	return 0;
}
