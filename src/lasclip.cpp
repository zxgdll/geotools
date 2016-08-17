/*
 * This tool clips a LiDAR point cloud to a vector shape.
 * Takes a shapefile and any number of LAS files as input.
 * Produces a single LAS file.
 *
 *  Author: Rob Skelly rob@dijital.ca
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <list>
#include <climits>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <ogr_api.h>
#include <ogrsf_frmts.h>
#include <gdal_priv.h>
#include <ogr_geometry.h>

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
#include "util.hpp"

namespace gg = geos::geom;

using namespace geotools::util;

namespace geotools {

	namespace las {

		void lasclip(std::string &outfile, std::string &shapefile, std::string &layername, 
			std::vector<std::string> &files, std::set<int> &classes, bool quiet) {

			if(outfile.empty()) 
				g_argerr("An output file is required.");
			if(shapefile.empty())
				g_argerr("A shape file is required.");
			if(files.size() == 0)
				g_argerr("At least one input file is required.");
			if(classes.size() == 0)
				g_warn("No classes specified, matching all classes.");

			/* Attempt to open and load geometries from the shape file. */
			OGRRegisterAll();
			OGRLayer *layer;
			OGRFeature *feat;
			OGRGeometry *og;
			OGRwkbGeometryType type;
			gg::GeometryCollection *geomColl;
			gg::Geometry *geom;

			OGRDataSource *ds = OGRSFDriverRegistrar::Open(shapefile.c_str(), FALSE);
			if(ds == nullptr)
				g_runerr("Couldn't open shapefile.");
			if(layername.empty()) {
				layer = ds->GetLayer(0);
			} else {
				layer = ds->GetLayerByName(layername.c_str());
			}
			if(layer == nullptr)
				g_runerr("Couldn't get layer.");

			type = layer->GetGeomType();
			if(type != wkbPolygon && type != wkbMultiPolygon)
				g_runerr("Geometry must be polygon or multipolygon.");

			const GEOSContextHandle_t gctx = OGRGeometry::createGEOSContext();
			const gg::GeometryFactory *gf = gg::GeometryFactory::getDefaultInstance();
			const gg::CoordinateSequenceFactory *cf = gf->getCoordinateSequenceFactory();
			std::vector<gg::Geometry *> geoms;

			while((feat = layer->GetNextFeature()) != NULL) {
				og = feat->GetGeometryRef();
				geom = (gg::Geometry *) og->exportToGEOS(gctx);
				geoms.push_back(geom);
			}

			GDALClose(ds);

			if(geoms.size() == 0)
				g_runerr("No geometries were found.");

			/* The geometry collection is used for checking whether a las file intersects
			   the region of interest. */
			geomColl = gf->createGeometryCollection(geoms);
			const gg::Envelope *env = geomColl->getEnvelopeInternal();
			Bounds cbounds(env->getMinX(), env->getMinY(), env->getMaxX(), env->getMaxY());

			/* Loop over files and figure out which ones are relevant. */
			liblas::ReaderFactory rf;
			liblas::Header *dsth = nullptr;
			std::vector<unsigned int> indices;

			for(unsigned int i = 0; i < files.size(); ++i) {

				std::ifstream in(files[i].c_str(), std::ios::in | std::ios::binary);
				liblas::Reader r = rf.CreateWithStream(in);
				liblas::Header h = r.GetHeader();

				if(i == 0)
					dsth = new liblas::Header(h);

				std::vector<gg::Coordinate> coords;
				coords.push_back(gg::Coordinate(h.GetMinX(), h.GetMinY()));
				coords.push_back(gg::Coordinate(h.GetMaxX(), h.GetMinY()));
				coords.push_back(gg::Coordinate(h.GetMaxX(), h.GetMaxY()));
				coords.push_back(gg::Coordinate(h.GetMinX(), h.GetMaxY()));
				coords.push_back(gg::Coordinate(h.GetMinX(), h.GetMinY()));
				
				gg::CoordinateSequence *cs = cf->create(&coords);
				gg::LinearRing *lr = gf->createLinearRing(cs);
				gg::Polygon *bounds = gf->createPolygon(lr, NULL);

				if(bounds->intersects(geomColl)) 
					indices.push_back(i);

				in.close();
			}

			if(indices.size() == 0) 
				g_runerr("No files matched the given bounds.");

			std::ofstream out(outfile, std::ios::out | std::ios::binary);
			liblas::WriterFactory wf;
			liblas::Writer w(out, *dsth);
			liblas::Header::RecordsByReturnArray recs;
			int count = 0;

			double bounds[] = { G_DBL_MAX_POS, G_DBL_MAX_NEG, G_DBL_MAX_POS, G_DBL_MAX_NEG, G_DBL_MAX_POS, G_DBL_MAX_NEG };

			g_trace("Using points from " << indices.size() << " files.");

			for(int i = 0; i < 5; ++i)
				recs.push_back(0);

		   	for(unsigned int i = 0; i < indices.size(); ++i) {

				std::ifstream in(files[indices[i]].c_str(), std::ios::in | std::ios::binary);
				liblas::Reader r = rf.CreateWithStream(in);
				liblas::Header h = r.GetHeader();

				g_trace("Processing file " << files[indices[i]]);

				while(r.ReadNextPoint()) {
					liblas::Point pt = r.GetPoint();
					
					int cls = pt.GetClassification().GetClass();
					if(classes.size() > 0 && !Util::inList(classes, cls)) 
						continue;
					
					double x = pt.GetX();
					double y = pt.GetY();
					const gg::Coordinate c(x, y);
					gg::Point *p = gf->createPoint(c);

					if(cbounds.contains(x, y) && geomColl->contains(p)) {
						++recs[cls];
						++count;
						w.WritePoint(pt);
						if(pt.GetX() < bounds[0]) bounds[0] = pt.GetX();
						if(pt.GetX() > bounds[1]) bounds[1] = pt.GetX();
						if(pt.GetY() < bounds[2]) bounds[2] = pt.GetY();
						if(pt.GetY() > bounds[3]) bounds[3] = pt.GetY();
						if(pt.GetZ() < bounds[4]) bounds[4] = pt.GetZ();
						if(pt.GetZ() > bounds[5]) bounds[5] = pt.GetZ();
					}
				}

				in.close();
			}

			// Set the total count and update the point record counts.
			dsth->SetMin(bounds[0], bounds[2], bounds[4]);
			dsth->SetMax(bounds[1], bounds[3], bounds[5]);
			dsth->SetPointRecordsCount(count);
			for(unsigned int i=0;i<recs.size();++i)
				dsth->SetPointRecordsByReturnCount(i, recs[i]);

			w.WriteHeader();

		}

	} // las

} // geotools

void usage() {
	std::cerr << "Usage: lasclip [options] srcfiles*\n"
	<< "  Produces a single LAS file containing points that are present in the shapefile.\n"
	<< " -o Speficy an output file.\n"
	<< " -s Specify a shapefile. Points inside any polygon will be preserved.\n"	
	<< " -l Specify a layer name.\n"	
	<< " -c A comma-delimited list of classes to accept.\n"	
	<< " -q Supress output.\n";
}

int main(int argc, char ** argv) {

	std::vector<std::string> files;
	std::string outfile;
	std::string shapefile;
	std::string layername;
	std::set<int> classes;
	bool quiet = false;

	/* Parse and check input. */

	for(int i=1;i<argc;++i) {
		std::string arg(argv[i]);
		if(arg == "-s") {
			shapefile = argv[++i];
		} else if(arg == "-o") {
			outfile = argv[++i];
		} else if(arg == "-l") {
			layername = argv[++i];
		} else if(arg == "-q") {
			quiet = true;
		} else if(arg == "-c") {
			Util::intSplit(classes, argv[++i]);
		} else {
			files.push_back(std::string(argv[i]));
		}
	}

	try {

		geotools::las::lasclip(outfile, shapefile, layername, files, classes, quiet);

	} catch(const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;

}
