/*
 * rastfit
 *
 * This program adjusts a raster's elevations to match another's using natural
 * neighbour interpolation or plane fitting. It does this by collecting a configurable
 * number of random samples constrained by an optional mask and using these as sites
 * to compute differences.
 *
 *  Created on: Jan 16, 2016
 *  Author: rob
 */

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <algorithm>

#include "util.hpp"
#include "Raster.hpp"
#include "interp/InterpPoint.hpp"
#include "interp/Interpolator.hpp"
#include "interp/IDWInterpolator.hpp"
#include "interp/AvgInterpolator.hpp"
#include "interp/PlanarInterpolator.hpp"
#include "interp/NaturalNeighbourInterpolator.hpp"
#include "interp/SimpleKrigingInterpolator.hpp"

//#include "ShapeWriter.hpp"

//ShapeWriter sw("/home/rob/Documents/geotools/data/out.shp");

using namespace interp;

double g_min(double a, double b) {
    return a > b ? b : a;
}

double g_max(double a, double b) {
    return a < b ? b : a;
}

double _random() {
    double r = ((double) std::rand()) / RAND_MAX;
    return r;
}

double g_sq(double a) {
    return a*a;
}

/**
 * Sorts points so that they're in row/col order, 
 * which improves the read efficiency on the raster.
 */
class samplesort {
public:

    samplesort(Raster<float> &r) : r(r) {
    }

    bool operator()(std::unique_ptr<InterpPoint> &a, std::unique_ptr<InterpPoint> &b) const {
        return (a.get(), b.get());
    }

    bool operator()(InterpPoint &a, InterpPoint &b) const {
        int ar = r.toBlockRow(a.y);
        int br = r.toBlockRow(b.y);
        int ac = r.toBlockCol(a.x);
        int bc = r.toBlockCol(b.x);
        if (ar == br) {
            return ac < bc;
        } else {
            return ar < br;
        }
    }
private:
    Raster<float> &r;
};

/**
 * Generates a set of sample points, limited to valid pixels in the mask.
 */
void generateMaskSamples(std::list<InterpPoint> &samples, Raster<char> &mask,
        Raster<float> &base, Raster<float> &adj, unsigned int numSamples) {
    std::vector<InterpPoint> pts;
    for (int r = 0; r < mask.rows(); ++r) {
        for (int c = 0; c < mask.cols(); ++c) {
            if (mask.get(c, r) == 1.0
                    && adj.isValid(mask.toX(c), mask.toY(r))
                    && base.isValid(mask.toX(c), mask.toY(r))) {
                InterpPoint p(mask.toX(c), mask.toY(r));
                p.z = adj.get(mask.toX(c), mask.toY(r)) - base.get(mask.toX(c), mask.toY(r));
                pts.push_back(p);
            }
        }
    }
    std::cerr << "Selecting from " << pts.size() << " points" << std::endl;
    std::random_shuffle(pts.begin(), pts.end());
    samples.assign(pts.begin(), pts.begin() + (numSamples > pts.size() ? pts.size() : numSamples));
}

/**
 * Generates a list of random samples within the bounds of a raster.
 */
void generateRandomSamples(std::list<InterpPoint> &samples, Raster<float> &adj, Raster<float> &base,
        double *bounds, unsigned int numSamples) {
    while (numSamples-- > 0) {
        double x = bounds[0] + (bounds[2] - bounds[0]) * _random();
        double y = bounds[1] + (bounds[3] - bounds[1]) * _random();
        if (adj.isValid(x, y) && base.isValid(x, y)) {
            InterpPoint p(x, y);
            p.z = adj.get(x, y) - base.get(x, y);
            samples.push_back(p);
        }
    }
}

/**
 * Print the sample locations in csv format with the original raster
 * difference, the adjustment value, the final value, and the residual.
 * Prints the variances at the bottom.
 */
void printSamples(std::list<InterpPoint> &samples) {
    std::cout << "x,y,base,adj,diff,final,resid" << std::endl;
    double d = 0.0, r = 0.0;
    for (auto it = samples.begin(); it != samples.end(); ++it) {
        d += g_sq(it->z);
        r += g_sq(it->resid);
        std::cout << std::setprecision(9)
                << it->x << "," << it->y << "," << it->z << "," << it->resid << std::endl;
    }
    std::cout << "Var before: " << (d / samples.size()) << ", after: " << (r / samples.size()) << std::endl;
}

/**
 * Compute the intersection of the bounds of the two rasters.
 */
void computeBounds(Raster<float> &a, Raster<float> &b, double *bounds) {
    bounds[0] = g_max(a.toX(0), b.toX(0));
    bounds[1] = g_max(a.toY(a.rows()), b.toY(b.rows()));
    bounds[2] = g_min(a.toX(a.cols()), b.toX(b.cols()));
    bounds[3] = g_min(a.toY(0), b.toY(0));
}

void adjust(std::string &basefile, std::string &adjfile,
        std::string &maskfile, std::string &outfile,
        interp::Interpolator &inter,
        unsigned int numSamples, double resolution, bool print) {

    if (basefile == outfile || adjfile == outfile || maskfile == outfile)
        throw "The output file must not be the same as any input file.";

    if (numSamples < 2)
        throw "Too few samples.";

    if (resolution <= 0.0)
        throw "Invalid resolution.";

    // Initializes source, destination rasters.
    Raster<float> base(basefile, 1, false);
    Raster<float> adj(adjfile, 1, false);
    std::string proj;
    adj.projection(proj);
    Raster<float> out(outfile, adj.minx(), adj.miny(), adj.maxx(), adj.maxy(), resolution, adj.nodata(), proj);

    base.get(0, 0);
    adj.get(0, 0);
    std::list<InterpPoint> samples;

    // Generate a list of samples, shuffle it and clip it to the
    // desired count.
    std::cerr << "Generating samples." << std::endl;
    if (maskfile.empty() || maskfile == "-") {
        double bounds[4];
        computeBounds(base, adj, bounds);
        generateRandomSamples(samples, adj, base, bounds, numSamples);
    } else {
        Raster<char> mask(maskfile, 1, false);
        generateMaskSamples(samples, mask, adj, base, numSamples);
    }

    if (samples.size() == 0)
        throw "No samples collected.";

    std::cerr << "Interpolating..." << std::endl;
    inter.interpolate(out, samples);

    if (print) {
        for (auto it = samples.begin(); it != samples.end(); ++it)
            it->resid = out.get(it->x, it->y);
        printSamples(samples);
    }

}

void usage() {
    std::cerr << "Usage: rastfit [options]" << std::endl
            << "This program will produce an adjustment raster to adjust one raster's " << std::endl
            << "elevations to match another's using a selected algorithm. The mask is " << std::endl
            << "optional and used to limit the placement of samples to where elevations are " << std::endl
            << "known to be good. The output is an adjustment raster with the same extent " << std::endl
            << "as the input raster, with the new resolution. " << std::endl
            << " -t  -- type            The type of adjustment. " << std::endl
            << "                        nn  - natural neighbours." << std::endl
            << "                        pl  - plane fit. " << std::endl
            << "                        avg - shift vertically by the average difference. " << std::endl
            << "                        idw - inverse distance weighting (use -e switch for exponent; default 1). " << std::endl
            << "                        sk  - Simple Kriging." << std::endl
            << " -m  -- mask            an optional raster whose non-zero pixels dictate the locations of allowable samples." << std::endl
            << " -r  -- resolution      the pixel size of the output." << std::endl
            << " -s  -- samples         the number of samples." << std::endl
            << " -b  -- base file       the raster to adjust the adjustment raster to." << std::endl
            << " -a  -- adjustment      file: the raster to adjust." << std::endl
            << " -o  -- output file     the output file." << std::endl
            << " -p  -- print           print a table with the points, samples, deltas and residuals." << std::endl
            << " -ie -- idw exponent    The IDW decay value." << std::endl
            << " -ip -- idw neighbours  Number of neighbours to consider. Leave out to consider all." << std::endl;
}

int main(int argc, char **argv) {

    try {

        std::string basefile;
        std::string adjfile;
        std::string maskfile;
        std::string outfile;
        std::string type;
        unsigned int numSamples = 0;
        double resolution = 0.0;
        bool print = false;
        double idwExp = 1.0;
        int idwNeigh = 0;

        for (int i = 0; i < argc; ++i) {
            std::string p(argv[i]);
            if (p == "-t") {
                type = argv[++i];
            } else if (p == "-m") {
                maskfile = argv[++i];
            } else if (p == "-s") {
                numSamples = atoi(argv[++i]);
            } else if (p == "-o") {
                outfile = argv[++i];
            } else if (p == "-a") {
                adjfile = argv[++i];
            } else if (p == "-r") {
                resolution = atof(argv[++i]);
            } else if (p == "-b") {
                basefile = argv[++i];
            } else if (p == "-p") {
                print = true;
            } else if (p == "-ie") {
                idwExp = atof(argv[++i]);
            } else if (p == "-ip") {
                idwNeigh = atoi(argv[++i]);
            }
        }

        if (type == "idw" && idwNeigh < 0)
            throw "IDW neighbour count must be >= 0.";
        if (numSamples <= 1)
            throw "Please try more than one sample.";
        if (resolution <= 0)
            throw "Invalid resolution.";
        if (basefile.empty())
            throw "Base file is required.";
        if (outfile.empty())
            throw "Output file is required.";
        if (adjfile.empty())
            throw "Adjustment file is required.";

        interp::Interpolator *inter;

        if (type == "nn") {
            inter = new interp::naturalneighbour::NaturalNeighbourInterpolator();
        } else if (type == "pl") {
            inter = new interp::planar::PlanarInterpolator();
        } else if (type == "avg") {
            inter = new interp::avg::AvgInterpolator();
        } else if (type == "idw") {
            inter = new interp::idw::IDWInterpolator(idwExp, idwNeigh);
        } else if (type == "sk") {
            interp::kriging::SimpleKrigingInterpolator *sk = new interp::kriging::SimpleKrigingInterpolator(argc, argv);
            inter = sk;
        } else {
            throw "Unknown interpolator type.";
        }

        adjust(basefile, adjfile, maskfile, outfile, *inter, numSamples, resolution, print);

        delete inter;

        std::cerr << "Done." << std::endl;

    } catch (const char *e) {
        std::cerr << e << std::endl;
        usage();
        return 1;
    }

    return 0;
}



