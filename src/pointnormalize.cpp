#include "pointnormalize.hpp"
#include "raster.hpp"
#include "pointstream.hpp"

using namespace geotools::point;
using namespace geotools::raster;
using namespace geotools::las;

void PointNormalize::normalize(const PointNormalizeConfig &config, const Callbacks *callbacks) {

	Raster<float> terrain(config.terrainFile, 1, false);
	Raster<float> chm(config.chmFile, 1, terrain);

	unsigned int fileIdx = 0;
	for(const std::string &pointFile : config.pointFiles) {

		if(callbacks)
			callbacks->overallCallback((fileIdx + 0.5f) / config.pointFiles.size());

		PointStream ps(pointFile);
		LASPoint pt;

		unsigned int ptIdx = 0;
		unsigned int ptInt = ps.pointCount() / 20;

		while(ps.next(pt)) {

			int col = terrain.toCol(pt.x);
			int row = terrain.toRow(pt.y);
			chm.set(col, row, pt.z - terrain.get(col, row));

			if(callbacks && ptIdx % ptInt == 0)
				callbacks->fileCallback((float) ptIdx / ps.pointCount());
			++ptIdx;
		}

		if(callbacks)
			callbacks->overallCallback((fileIdx + 1.0f) / config.pointFiles.size());
		++fileIdx;
	}	

	if(callbacks) {
		callbacks->overallCallback(1.0f);
		callbacks->fileCallback(1.0f);
	}
}