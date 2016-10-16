#include "liblas/liblas.hpp"
#include "boost/filesystem.hpp"

#include "pointnormalize.hpp"
#include "raster.hpp"
#include "pointstream.hpp"
#include "util.hpp"

using namespace geotools::point;
using namespace geotools::raster;
using namespace geotools::las;

void PointNormalize::normalize(const PointNormalizeConfig &config, const Callbacks *callbacks) {

	if(config.pointOutputDir.empty())
		g_argerr("A point output dir must be provided.");

 	if(config.threads > 0) {
 		g_debug(" -- pointstats running with " << config.threads << " threads");
 		omp_set_dynamic(1);
 		omp_set_num_threads(config.threads);
 	} else {
 		g_argerr("Run with >=1 threads.");
 	}	

	if(!Util::mkdir(config.pointOutputDir))
		g_argerr("Couldn't create output dir " << config.pointOutputDir);

 	using namespace boost::filesystem;

	path outDir(config.pointOutputDir);
	liblas::ReaderFactory rf;

	std::vector<std::string> files(config.pointFiles.begin(), config.pointFiles.end());

	bool overwrite = true; // TODO: Configurable!

	#pragma omp parallel
	{

		Raster<float> terrain(config.terrainFile, 1, false);

		#pragma omp for
		for(unsigned int i = 0; i < files.size(); ++i) {
			
			const std::string &pointFile = files[i];

			path oldPath(pointFile);
			path newPath(outDir / oldPath.filename());
			
			g_debug(" -- normalize (point) processing " << pointFile << " to " << newPath);

			if(!overwrite && exists(newPath)) {
				g_warn("The output file " << newPath << " already exists.");
				continue;
			}

			std::ifstream instr(pointFile.c_str(), std::ios::binary);
			liblas::Reader lasReader = rf.CreateWithStream(instr);
			liblas::Header lasHeader = lasReader.GetHeader();

			std::ofstream outstr(newPath.c_str(), std::ios::binary);
			liblas::Writer lasWriter(outstr, lasHeader);

			while(lasReader.ReadNextPoint()) {
				const liblas::Point &opt = lasReader.GetPoint();
				int col = terrain.toCol(opt.GetX());
				int row = terrain.toRow(opt.GetY());
				if(terrain.has(col, row)) {
					liblas::Point npt(opt);
					if(!terrain.isNoData(col, row)) {
						npt.SetZ(opt.GetZ() - terrain.get(col, row));
						lasWriter.WritePoint(npt);
					}
				}
			}
			
			instr.close();
			outstr.close();
		}

	}

}