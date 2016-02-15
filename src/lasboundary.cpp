/*
 * lasboundary.cpp
 *
 * This tool creates a shapefile representing the boundary of a LiDAR
 * point cloud as represented by an alpha shape. See CGAL
 * alpha shape docs for more info.
 *
 *  Created on: Mar 13, 2015
 *      Author: rob
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <list>
#include <climits>
#include <memory>
#include <cstring>
#include <math.h>

#include <omp.h>
#include <ogr_spatialref.h>
#include <ogrsf_frmts.h>
#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <liblas/liblas.hpp>

#include "Util.hpp"

#define LAS_EXT ".las"

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;
namespace las = liblas;

/**
 * Gets a list of the files in the directory pointed to by srcDir, 
 * or returns a list containing srcDir if srcDir is a file.
 */
void listFiles(std::vector<std::string> &files, const char *srcDir) {
	fs::path srcdir_p(srcDir);
	if(!exists(srcdir_p))
		throw std::invalid_argument("Path not found.");
	if(is_directory(srcdir_p)) {
		fs::directory_iterator end;
		for(fs::directory_iterator iter(srcdir_p); iter != end; iter++) {
			std::string fname = iter->path().string();
			std::string lname(fname);
			alg::to_lower(lname);
			if(alg::ends_with(lname, LAS_EXT)) 
				files.push_back(fname);
		}
	} else {
		std::string s(srcDir);
		files.push_back(s);
	}
}

/**
 * Splits a string of comma-delimited integers into a set.
 */
void intSplit(std::set<int> &classes, const char *classStr) {
	std::stringstream ss(classStr);
	std::string item;
	while(std::getline(ss, item, ','))
		classes.insert(atoi(item.c_str()));
}

/**
 * Returns true if the given int is in the set of ints.
 */
bool inList(std::set<int> &classes, int cls) {
	return classes.find(cls) != classes.end();
}

/**
 * Returns the square of a given float.
 */
inline float _sq(float a) {
	return a * a;
}

void usage() {
	std::cerr << "Usage: lasboundary [options] -i <src dir> -o <dst file>" << std::endl;
	std::cerr << "	This program creates a Shapefile containing the boundary " << std::endl;
	std::cerr << "  of a point cloud represented by a set of LAS files. The " << std::endl;
	std::cerr << "  boundary is an alpha shape based on the Delaunay triangulation  " << std::endl;
	std::cerr << "  with alpha as the square of the given radius." << std::endl;
	std::cerr << "  src dir  - The source directory or a single LAS file." << std::endl;
	std::cerr << "  dst file - The output Shapefile." << std::endl;
	std::cerr << "  -c       - A comma-delimited list of integers indicating " << std::endl;
	std::cerr << "             which classes to preserve (e.g. 2 = ground). Defaults to all." << std::endl;
	std::cerr << "  -r       - The resolution of the grid. Default 1.0m." << std::endl;
	std::cerr << "  -s       - The integer EPSG ID of the coordinate reference system." << std::endl;
}

int main(int argc, char **argv) {

	if(argc < 6) {
		usage();
		return 1;
	}

	char *srcDir = 0;
	char *dstFile = 0;
	int srid = 0;
	double res = 1.0;
	std::set<int> classes;
	bool hasClasses = false;

	for(int i = 0; i < argc; ++i) {
		std::string s(argv[i]);
		if(s == "-c") {
			// Gets the set of classes to keep
			intSplit(classes, argv[++i]);
			hasClasses = true;
		} else if(s == "-r") {
			 res = atof(argv[++i]);
		} else if(s == "-s") {
			srid = atoi(argv[++i]);
		} else if(s == "-i") {
			srcDir = argv[++i];
		} else if(s == "-o") {
			dstFile = argv[++i];
		}
	}

	if(!srcDir) {
		std::cerr << "No source dir given." << std::endl;
		usage();
		return 1;
	}

	if(!dstFile) {
		std::cerr << "No dest file given." << std::endl;
		usage();
		return 1;
	}

	if(srid <= 0) {
		std::cerr << "No SRID given." << std::endl;
		usage();
		return 1;
	}

	// Gets the list of files from the source dir.
	std::vector<std::string> files;
	listFiles(files, srcDir);

	if(files.size() == 0) {
		std::cerr << "No files found." << std::endl;
		usage();
		return 1;
	}

	las::ReaderFactory rf;
	double bounds[4] = { FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX } ;
	double bounds0[4] = { FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX } ;

	for(unsigned int i=0; i<files.size(); ++i) {
		// Open the file and create a reader.
		std::cerr << "Checking " << files[i] << std::endl;
		std::ifstream in(files[i].c_str());
		las::Reader reader = rf.CreateWithStream(in);
		las::Header header = reader.GetHeader();
		Util::computeLasBounds(header, bounds0, 2);
		Util::expand(bounds, bounds0, 2);
		in.close();
	}

	Util::snapBounds(bounds, res, 2);
	int cols = (int) ((bounds[2] - bounds[0]) + res);
	int rows = (int) ((bounds[3] - bounds[1]) + res);
	double width = bounds[2] - bounds[0];
	double height = bounds[3] - bounds[1];
	int col, row;
	char *grid = (char *) calloc(cols * rows, sizeof(char));

	for(unsigned int i=0; i<files.size(); ++i) {
		// Open the file and create a reader.
		std::cerr << "Processing " << files[i] << std::endl;
		std::ifstream in(files[i].c_str());
		las::Reader reader = rf.CreateWithStream(in);
		las::Header header = reader.GetHeader();
		while(reader.ReadNextPoint()) {
			las::Point pt = reader.GetPoint();
			if(!hasClasses || inList(classes, pt.GetClassification().GetClass())) {
				col = (int) ((pt.GetX() - bounds[0]) / width * res * cols);
				row = (int) ((pt.GetY() - bounds[1]) / height * res * rows);
				grid[row * cols + col] = 1;
			}
		}
	}

	// Build the shapefile output.
	GDALAllRegister();
	OGRRegisterAll();
	double transform[6] = {bounds[0], res, 0, bounds[3], 0, -res};
	const char *shpFormat = "ESRI Shapefile";
	const char *tifFormat = "GTiff";

    OGRSpatialReference ref;
    ref.importFromEPSG(srid);
    char *crs = (char *) malloc(1024 * sizeof(char));
    ref.exportToWkt(&crs);

    OGRSFDriver *shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(shpFormat);
    GDALDriver *tifDriver = GetGDALDriverManager()->GetDriverByName(tifFormat);

    GDALDataset *tifDS = tifDriver->Create("/tmp/lasboundary.tif", cols, rows, 1, GDT_Byte, 0);
    tifDS->SetGeoTransform(transform);
    tifDS->SetProjection(crs);
    tifDS->RasterIO(GF_Write, 0, 0, cols, rows, grid, cols, rows, GDT_Byte, 1, 0, 0, 0, 0);
    OGRFree(crs);

    OGRDataSource *shpDS = shpDriver->CreateDataSource(dstFile, NULL);
    OGRLayer *layer = shpDS->CreateLayer("boundary", &ref, wkbMultiPolygon, NULL);
	OGRFieldDefn fvalue("value", OFTReal);
	if(layer->CreateField(&fvalue) != OGRERR_NONE) {
		std::cerr << "Failed to create value field." << std::endl;
		return 1;
	}

    //GDALPolygonize(tifDS, NULL, layer, 0, NULL, NULL, NULL);
	std::cerr << "Writing to file." << std::endl;
    GDALClose(tifDS);
    GDALClose(shpDS);

    return 0;

}
