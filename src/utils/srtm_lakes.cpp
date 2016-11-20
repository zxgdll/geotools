#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <queue>
#include <set>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include "Raster.hpp"
#include "util.hpp"

unsigned long _join(int a, int b) {
    return (unsigned long) a << 32 | b;
}

void _split(unsigned long j, int *a, int *b) {
    *a = (int) (j >> 32);
    *b = (int) (j & 0xffffffff);
}

void delineate(int col, int row, Grid<float> &grid, Grid<float> &water, float nodata, unsigned int minsize) {
    float el = grid(col, row);
    if (el == nodata)
        return;
    std::queue<unsigned long> cells;
    std::set<unsigned long> out;
    int c, r;
    cells.push(_join(col, row));
    while (cells.size() > 0) {
        unsigned long i = cells.front();
        cells.pop();
        _split(i, &c, &r);
        if (c < 0 || r < 0 || c >= grid.cols() || r >= grid.rows())
            continue;
        float v = grid(c, r);
        if (v != nodata) {
            grid(c, r, nodata);
            if (v == el) {
                out.insert(i);
                cells.push(_join(c - 1, r));
                cells.push(_join(c + 1, r));
                cells.push(_join(c, r - 1));
                cells.push(_join(c, r + 1));
                cells.push(_join(c - 1, r - 1));
                cells.push(_join(c - 1, r + 1));
                cells.push(_join(c + 1, r - 1));
                cells.push(_join(c + 1, r + 1));
            }
        }
    }
    if (out.size() >= minsize) {
        int c, r;
        for (auto it = out.begin(); it != out.end(); ++it) {
            _split(*it, &c, &r);
            water(c, r, el);
        }
    }
}

int main(int argc, char **argv) {

    try {

        if (argc < 4)
            throw "Usage: strm_lakes <infile> <outfile> <minsize>";

        std::string infile = argv[1];
        std::string outfile = argv[2];
        unsigned int minsize = atoi(argv[3]);

        Raster<float> inrast(infile);
        std::string proj;
        inrast.projection(proj);
        double nodata = inrast.nodata();

        std::cerr << "Raster: " << inrast.cols() << ", " << inrast.rows() << std::endl;

        Grid<float> grid(inrast.cols(), inrast.rows());
        grid.fill(nodata);

        Grid<float> water(inrast.cols(), inrast.rows());

        inrast.readBlock(grid);

        for (int r = 0; r < inrast.rows(); ++r) {
            std::cerr << "Row " << r << " of " << inrast.rows() << std::endl;
            for (int c = 0; c < inrast.cols(); ++c) {
                delineate(c, r, grid, water, nodata, minsize);
            }
        }

        Raster<float> outrast(outfile, inrast);
        outrast.writeBlock(water);

    } catch (const char *e) {
        std::cerr << e << std::endl;
        return 1;
    }
    return 0;
}