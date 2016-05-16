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


#include <geos/triangulate/DelaunayTriangulationBuilder.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Point.h>
#include <geos/geom/Coordinate.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <liblas/liblas.hpp>

#include "geotools.h"
#include "Util.hpp"
#include "Raster.hpp"
#include "Vector.hpp"

#define LAS_EXT ".las"

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;
namespace las = liblas;

/**
 * Gets a list of the files in the directory pointed to by srcDir, 
 * or returns a list containing srcDir if srcDir is a file.
 */
void listFiles(std::vector<std::string> &files, std::string &srcDir) {
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

bool fullNeighbours(Grid<char> &grid, int col, int row) {
	if(col == 0 || row == 0 || col >= grid.cols() - 1 || row >= grid.rows() - 1)
		return false;
	if(
		!grid.get(col-1, row-1) || 
		!grid.get(col, row-1) || 
		!grid.get(col+1, row-1) ||
		!grid.get(col-1, row) || 
		!grid.get(col+1, row) ||
		!grid.get(col-1, row+1) || 
		!grid.get(col, row+1) || 
		!grid.get(col+1, row+1))
		return false;
	return true;
}

void buildBoundary(std::string &srcDir, std::string &dstFile, int srid, double res, std::set<int> &classes) {

	if(srcDir.empty()) 
		throw "No source dir given.";

	if(dstFile.empty())
		throw "No dest file given.";

	if(srid <= 0)
		throw "No SRID given.";

	// Gets the list of files from the source dir.
	std::vector<std::string> files;
	listFiles(files, srcDir);

	if(files.size() == 0)
		throw "No files found.";

	las::ReaderFactory rf;
	std::vector<double> bounds = { DBL_MAX_POS, DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_NEG } ;

	for(unsigned int i=0; i<files.size(); ++i) {
		// Open the file and create a reader.
		std::cerr << "Checking " << files[i] << std::endl;
		std::ifstream in(files[i].c_str());
		las::Reader reader = rf.CreateWithStream(in);
		las::Header header = reader.GetHeader();
		std::vector<double> bounds0 = { DBL_MAX_POS, DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_NEG } ;
		Util::computeLasBounds(header, bounds0, 2);
		Util::expand(bounds, bounds0, 2);
		in.close();
	}

	Util::snapBounds(bounds, res, 2);
	int cols = (int) ((bounds[2] - bounds[0]) + res);
	int rows = (int) ((bounds[3] - bounds[1]) + res);
	double width = bounds[2] - bounds[0];
	double height = bounds[3] - bounds[1];

	MemRaster<char> grid(cols, rows);

	for(unsigned int i=0; i<files.size(); ++i) {
		// Open the file and create a reader.
		std::cerr << "Processing " << files[i] << std::endl;
		std::ifstream in(files[i].c_str());
		las::Reader reader = rf.CreateWithStream(in);
		las::Header header = reader.GetHeader();
		while(reader.ReadNextPoint()) {
			las::Point pt = reader.GetPoint();
			if(classes.size() == 0 || inList(classes, pt.GetClassification().GetClass())) {
				int col = (int) ((pt.GetX() - bounds[0]) / res);
				int row = (int) ((pt.GetY() - bounds[1]) / res);
				grid.set(col, row, 1);
			}
		}
	}

	using namespace geos::geom;
	using namespace geos::triangulate;

	GeometryFactory::unique_ptr gf = GeometryFactory::create();

	std::vector<Geometry*> coords;

	// Create a point set from the grid cells with returns in them.
	for(int row = 0; row < grid.rows(); ++row) {
		for(int col = 0; col < grid.cols(); ++col) {
			if(grid.get(col, row) && !fullNeighbours(grid, col, row))
				coords.push_back(gf->createPoint(Coordinate(col * res + bounds[0] + res / 2.0, row * res + bounds[1] - res / 2.0)));
		}
	}

	// Build Delaunay.
	GeometryCollection *mp = gf->createGeometryCollection(coords);
	DelaunayTriangulationBuilder dtb;
	dtb.setSites(*mp);

	// Get the edges.
	std::unique_ptr<MultiLineString> boundary = dtb.getEdges(*gf);

	MultiLineString *edges = boundary.get();
	std::vector<Geometry*> newlines;
	for(unsigned int i = 0; i < edges->getNumGeometries(); ++i) {
		Geometry *line = (Geometry *) edges->getGeometryN(i);
		if(line->getLength() < 10.0) {
			newlines.push_back(line);
		}
	}

	MultiLineString *newbounds = gf->createMultiLineString(&newlines);

	// Build the vector.
	std::map<std::string, int> attribs;
	attribs["value"] = Vector::INTEGER;
	Vector vec(dstFile, Vector::MULTILINE, std::string("epsg:26912"), attribs);

	std::unique_ptr<Geom> g = vec.addMultiLine(*newbounds);
	g->setAttribute("value", 1);

}

int main(int argc, char **argv) {

	if(argc < 6) {
		usage();
		return 1;
	}

	std::string srcDir;
	std::string dstFile;
	int srid = 0;
	double res = 1.0;
	std::set<int> classes;

	for(int i = 0; i < argc; ++i) {
		std::string s(argv[i]);
		if(s == "-c") {
			// Gets the set of classes to keep
			intSplit(classes, argv[++i]);
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

	try {

		buildBoundary(srcDir, dstFile, srid, res, classes);
	} catch(const char *e) {
		std::cerr << e << std::endl;
		usage();
		return 1;
	}

    return 0;

}
