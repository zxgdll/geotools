/*
 *  Created on: Feb 22, 2016
 *  Author: rob
 */

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <algorithm>

#include "geotools.h"
#include "util.hpp"
#include "Raster.hpp"
#include "interp/InterpPoint.hpp"
#include "interp/Interpolator.hpp"
#include "interp/IDWInterpolator.hpp"
#include "interp/AvgInterpolator.hpp"
#include "interp/PlanarInterpolator.hpp"
#include "interp/NaturalNeighbourInterpolator.hpp"
#include "interp/SimpleKrigingInterpolator.hpp"

double _random() {
    return ((double) std::rand()) / RAND_MAX;
}

/**
 * Load the samples from a csv file. The file must have x, y and z headers.
 */
void loadSamples(std::string &datafile, std::list<interp::InterpPoint> &samples) {
    g_runerr("CSV read not implemented.")
}

/**
 * Interpolate points from the data file into outfile, using templatefile as
 * the template. The parameterized interpolator performs the interpolation. The
 * resolution overrides the template's resolution.
 */
void interpolate(std::string &datafile, std::string &templatefile, std::string &outfile,
        interp::Interpolator &inter,
        double resolution, bool print) {

    if (templatefile == outfile)
        throw "The output file must not be the same as the template file.";

    if (resolution <= 0.0)
        throw "Invalid resolution.";

    // Initializes source, destination rasters.
    Raster<float> tpl(templatefile, 1, false);
    std::string proj;
    tpl.projection(proj);
    Raster<float> out(outfile, tpl.minx(), tpl.miny(), tpl.maxx(), tpl.maxy(),
            resolution, tpl.nodata(), proj);

    std::list<interp::InterpPoint> samples;
    loadSamples(datafile, samples);

    if (samples.size() == 0)
        throw "No samples loaded.";

    // for(auto it = samples.begin(); it != samples.end(); ++it)
    //	std::cerr << "sample " << it->x << ", " << it->y << std::endl;

    std::cerr << "Interpolating..." << std::endl;
    inter.interpolate(out, samples);

    if (print) {
        std::cerr << "You asked for it." << std::endl;
    }
}

void usage() {
    std::cerr << "Usage: interp [options]" << std::endl
            << " -t  -- type            The type of adjustment. " << std::endl
            << "                         nn  - natural neighbours." << std::endl
            << "                         pl  - plane fit. " << std::endl
            << "                         avg - shift vertically by the average difference. " << std::endl
            << "                         idw - inverse distance weighting (use -e switch for exponent; default 1). " << std::endl
            << "                         sk  - Simple Kriging." << std::endl
            << "                         lo  - LOESS Smoothing." << std::endl
            << " -r  -- resolution      The pixel size of the output." << std::endl
            << " -i  -- template file   A template file to produce the output file." << std::endl
            << " -o  -- output file     The output file." << std::endl
            << " -d  -- data file       A CSV file with data." << std::endl
            << " -ie -- idw exponent    The IDW decay value." << std::endl
            << " -ip -- idw neighbours  Number of neighbours to consider. Leave out to consider all." << std::endl
            << " -an -- avg neighbours  Number of neighbours for averaging. Default all." << std::endl;
}

int main(int argc, char **argv) {

    try {

        std::string templatefile;
        std::string outfile;
        std::string datafile;
        std::string type;
        double resolution = 0.0;
        bool print = false;
        double idwExp = 1.0;
        int idwNeigh = 0;
        unsigned int avgNeigh = 0;

        for (int i = 0; i < argc; ++i) {
            std::string p(argv[i]);
            if (p == "-t") {
                type = argv[++i];
            } else if (p == "-o") {
                outfile = argv[++i];
            } else if (p == "-i") {
                templatefile = argv[++i];
            } else if (p == "-r") {
                resolution = atof(argv[++i]);
            } else if (p == "-p") {
                print = true;
            } else if (p == "-ie") {
                idwExp = atof(argv[++i]);
            } else if (p == "-ip") {
                idwNeigh = atoi(argv[++i]);
            } else if (p == "-d") {
                datafile = argv[++i];
            } else if (p == "-an") {
                avgNeigh = atoi(argv[++i]);
            }
        }

        if (type == "idw" && idwNeigh < 0)
            throw "IDW neighbour count must be >= 0.";
        if (resolution <= 0)
            throw "Invalid resolution.";
        if (outfile.empty())
            throw "An output file is required.";
        if (templatefile.empty())
            throw "A template file is required.";
        if (datafile.empty())
            throw "A data file is required.";

        interp::Interpolator *inter;

        if (type == "nn") {
            inter = new interp::naturalneighbour::NaturalNeighbourInterpolator();
        } else if (type == "pl") {
            inter = new interp::planar::PlanarInterpolator();
        } else if (type == "avg") {
            inter = new interp::avg::AvgInterpolator(avgNeigh);
        } else if (type == "idw") {
            inter = new interp::idw::IDWInterpolator(idwExp, idwNeigh);
        } else if (type == "sk") {
            inter = new interp::kriging::SimpleKrigingInterpolator(argc, argv);
        } else {
            throw "Unknown interpolator type.";
        }

        interpolate(datafile, templatefile, outfile, *inter, resolution, print);

        delete inter;

        std::cerr << "Done." << std::endl;

    } catch (const char *e) {
        std::cerr << e << std::endl;
        usage();
        return 1;
    }

    return 0;
}



