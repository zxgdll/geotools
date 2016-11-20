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
#include "util.hpp"
#include "lasutil.hpp"
#include "raster.hpp"
#include "Vector.hpp"

#define LAS_EXT ".las"

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;
namespace geom = geos::geom;


#include "raster.cpp"


using namespace geotools::util;
using namespace geotools::raster;
using namespace geotools::vector;

namespace geotools {

    namespace las {

        namespace util {

            /**
             * Gets a list of the files in the directory pointed to by srcDir, 
             * or returns a list containing srcDir if srcDir is a file.
             */
            void listFiles(std::vector<std::string> &files, std::string &srcDir) {
                fs::path srcdir_p(srcDir);
                if (!exists(srcdir_p))
                    g_argerr("Path not found.");
                if (is_directory(srcdir_p)) {
                    fs::directory_iterator end;
                    for (fs::directory_iterator iter(srcdir_p); iter != end; iter++) {
                        std::string fname = iter->path().string();
                        std::string lname(fname);
                        alg::to_lower(lname);
                        if (alg::ends_with(lname, LAS_EXT))
                            files.push_back(fname);
                    }
                } else {
                    std::string s(srcDir);
                    files.push_back(s);
                }
            }

            bool fullNeighbours(Grid<char> &grid, int col, int row) {
                if (col == 0 || row == 0 || col >= grid.cols() - 1 || row >= grid.rows() - 1)
                    return false;
                if (
                        !grid.get(col - 1, row - 1) ||
                        !grid.get(col, row - 1) ||
                        !grid.get(col + 1, row - 1) ||
                        !grid.get(col - 1, row) ||
                        !grid.get(col + 1, row) ||
                        !grid.get(col - 1, row + 1) ||
                        !grid.get(col, row + 1) ||
                        !grid.get(col + 1, row + 1))
                    return false;
                return true;
            }

        } // itil

        void lasboundary(std::string &srcDir, std::string &dstFile, int srid, double res, std::set<int> &classes) {

            if (srcDir.empty())
                g_argerr("No source dir given.");
            if (dstFile.empty())
                g_argerr("No dest file given.");
            if (srid <= 0)
                g_argerr("No SRID given.");

            using namespace geotools::las::util;

            // Gets the list of files from the source dir.
            std::vector<std::string> files;
            listFiles(files, srcDir);

            if (files.size() == 0)
                g_argerr("No files found.");

            liblas::ReaderFactory rf;
            Bounds bounds;
            bounds.collapse();

            for (unsigned int i = 0; i < files.size(); ++i) {
                // Open the file and create a reader.
                g_trace("Checking " << files[i]);
                std::ifstream in(files[i].c_str());
                liblas::Reader reader = rf.CreateWithStream(in);
                liblas::Header header = reader.GetHeader();
                Bounds bounds0;
                bounds0.collapse();
                if (!LasUtil::computeLasBounds(header, bounds0, 2))
                    LasUtil::computeLasBounds(reader, bounds0, 2);
                bounds.extend(bounds0);
                in.close();
            }

            bounds.snap(res);
            int cols = (int) ((bounds[2] - bounds[0]) + res);
            int rows = (int) ((bounds[3] - bounds[1]) + res);
            double width = bounds[2] - bounds[0];
            double height = bounds[3] - bounds[1];

            MemRaster<char> grid(cols, rows);

            for (unsigned int i = 0; i < files.size(); ++i) {
                // Open the file and create a reader.
                g_trace("Processing " << files[i]);
                std::ifstream in(files[i].c_str());
                liblas::Reader reader = rf.CreateWithStream(in);
                liblas::Header header = reader.GetHeader();
                while (reader.ReadNextPoint()) {
                    liblas::Point pt = reader.GetPoint();
                    if (classes.size() == 0 || Util::inList(classes, pt.GetClassification().GetClass())) {
                        int col = (int) ((pt.GetX() - bounds[0]) / res);
                        int row = (int) ((pt.GetY() - bounds[1]) / res);
                        grid.set(col, row, 1);
                    }
                }
            }

            using namespace geos::geom;
            using namespace geos::triangulate;

            const geom::GeometryFactory *gf = geom::GeometryFactory::getDefaultInstance();
            std::vector<Geometry*> coords;

            // Create a point set from the grid cells with returns in them.
            for (int row = 0; row < grid.rows(); ++row) {
                for (int col = 0; col < grid.cols(); ++col) {
                    if (grid.get(col, row) && !fullNeighbours(grid, col, row))
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
            for (unsigned int i = 0; i < edges->getNumGeometries(); ++i) {
                Geometry *line = (Geometry *) edges->getGeometryN(i);
                if (line->getLength() < 10.0) {
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

    } // las

} // geotools

void usage() {
    std::cerr << "Usage: lasboundary [options] -i <src dir> -o <dst file>\n"
            << "	This program creates a Shapefile containing the boundary \n"
            << "  of a point cloud represented by a set of LAS files. The \n"
            << "  boundary is an alpha shape based on the Delaunay triangulation  \n"
            << "  with alpha as the square of the given radius.\n"
            << "  src dir  - The source directory or a single LAS file.\n"
            << "  dst file - The output Shapefile.\n"
            << "  -c       - A comma-delimited list of integers indicating \n"
            << "             which classes to preserve (e.g. 2 = ground). Defaults to all.\n"
            << "  -r       - The resolution of the grid. Default 1.0m.\n"
            << "  -s       - The integer EPSG ID of the coordinate reference system.\n";
}

int main(int argc, char **argv) {

    std::string srcDir;
    std::string dstFile;
    int srid = 0;
    double res = 1.0;
    std::set<int> classes;

    for (int i = 0; i < argc; ++i) {
        std::string s(argv[i]);
        if (s == "-c") {
            // Gets the set of classes to keep
            Util::intSplit(classes, argv[++i]);
        } else if (s == "-r") {
            res = atof(argv[++i]);
        } else if (s == "-s") {
            srid = atoi(argv[++i]);
        } else if (s == "-i") {
            srcDir = argv[++i];
        } else if (s == "-o") {
            dstFile = argv[++i];
        }
    }

    try {

        geotools::las::lasboundary(srcDir, dstFile, srid, res, classes);

    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        usage();
        return 1;
    }

    return 0;

}
