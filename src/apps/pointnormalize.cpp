#include "pointnormalize.hpp"

void usage() {
	std::cerr << "Usage: pointnormalize <terrain file> <chmfile> <point file [point file [point file ...]]>\n"
		<< " -v                          Verbose output.\n"
		<< " -h                          Print this message.\n"
		<< " --threads                   The number of threads to use for computing output.\n"
		<< " -gui                        Run the graphical user interface.\n";
}

int main(int argc, char **argv) {
	
int main(int argc, char **argv) {

	try {

		std::string terrainFile;
		std::string chmFile;
		std::list<std::string> pointFiles;
		int threads = 1;

		g_loglevel(0);
		
		for(int i = 1; i < argc; ++i) {
			std::string s(argv[i]);
			if(s == "-h") {
				usage();
				return 0;
			} else if(s == "-gui") {
				gui = true;
			} else if(s == "-v") {
				g_loglevel(G_LOG_DEBUG);
			} else if(s == "--threads") {
				threads = atoi(argv[++i]);
			} else {
				if(terrainFile.empty()) {
					terrainFile.assign(s);
				} else if(chmFile.empty()) {
					chmFile.assign(s);
				} else {
					pointFiles.push_back(argv[i]);
				}
			}
		}

		if(gui) {
			return runWithUI(argc, argv);
		} else {
			PointNormaize pn;
			PointNormailzeConfig config;
			config.terrainFile = terrainFile;
			config.chmFile = chmFile;
			config.pointFiles = pointFiles;
			config.threads = threads;
			
			pn.pointnormalize(config);
		}

	} catch(const std::exception &ex) {
		g_error(ex.what());
		usage();
		return 1;
	}

	return 0;
}
