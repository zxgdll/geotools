/*
 * Grids a point cloud represented by one or more LAS files.
 * Can produce grids of from intensity and elevation, using
 * minimum, maximum, mean, std dev, density, variance and count.
 *
 * Authored by: Rob Skelly rob@dijital.ca
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <list>
#include <climits>
#include <memory>
#include <cstring>
#include <math.h>
#include <exception>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <liblas/liblas.hpp>

#include "geotools.h"
#include "Util.hpp"
#include "Raster.hpp"

#define TYPE_MIN 1
#define TYPE_MAX 2
#define TYPE_MEAN 3
#define TYPE_DENSITY 4
#define TYPE_VARIANCE 7
#define TYPE_STDDEV 8
#define TYPE_PVARIANCE 12
#define TYPE_PSTDDEV 13
#define TYPE_COUNT 9
#define TYPE_QUANTILE 10
#define TYPE_MEDIAN 11

#define LAS_EXT ".las"

#define ATT_HEIGHT 1
#define ATT_INTENSITY 2

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;
namespace las = liblas;

/**
 * Interpret the value of a string attribute name, return the constant int value.
 */
int parseAtt(char *attStr) {
	if(!strcmp("intensity", attStr)) {
		return ATT_INTENSITY;
	} else if(!strcmp("height", attStr)) {
		return ATT_HEIGHT;
	} 
	return 0;
}

/**
 * Interpret the output type and return the constant int value.
 */
int parseType(char *typeStr) {
	if(!strcmp("min", typeStr)) {
		return TYPE_MIN;
	} else if(!strcmp("max", typeStr)) {
		return TYPE_MAX;
	} else if(!strcmp("mean", typeStr)) {
		return TYPE_MEAN;
	} else if(!strcmp("density", typeStr)) {
		return TYPE_DENSITY;
	} else if(!strcmp("variance", typeStr)) {
		return TYPE_VARIANCE;
	} else if(!strcmp("stddev", typeStr)) {
		return TYPE_STDDEV;
	} else if(!strcmp("pvariance", typeStr)) {
		return TYPE_PVARIANCE;
	} else if(!strcmp("pstddev", typeStr)) {
		return TYPE_PSTDDEV;
	} else if(!strcmp("count", typeStr)) {
		return TYPE_COUNT;
	} else if(!strcmp("median", typeStr)) {
		return TYPE_MEDIAN;
	}
	return 0;
}

/**
 * Comparator for sorting doubles.
 */
int _fcmp(const void * a, const void * b) {
	const double * aa = (const double *) a;
	const double * bb = (const double *) b;
	return (*aa > *bb) - (*bb > *aa);
}

void usage() {
	_print("Usage: lasgrid <options> <file [file [file]]>\n"
				<< " -o <output file>\n"
				<< " -t <type>                   Output median, mean, max, min, variance (sample), pvariance (population),\n"
			 	<< "                             count, density, stddev (sample), pstddev (population). Default mean.\n"
				<< " -r <resolution>             Resolution (default 2).\n"
				<< " -s <srid>                   The EPSG ID of the CRS.\n"
				<< " -c <classes>                Comma-delimited (e.g. '2,0' (ground and unclassified)).\n"
				<< " -a <attribute>              Use height, intensity (default height).\n"
				<< " -d <radius>                 Radius (not diameter); use zero for cell bounds.\n"
				<< "                             For example, if the cell size is 2, the circumcircle's radius is sqrt(2) (~1.41).\n"
				<< " -b <minx miny maxx maxy>    Extract points from the given box and create a raster of this size.\n"
				<< " -f                          Fill voids.\n"
				<< " -v                          Verbose output.");
}

void vector_dealloc(std::vector<double> *item) {
	delete item;
}

/**
 * Returns true if the point is within the radius associated with a cell's centroid.
 * @param px The x coordinate of the point.
 * @param py The y coordinate of the point.
 * @param col The column of the cell of interest.
 * @param row The row of the cell of interest.
 * @param radius The radius around the cell's centroid.
 * @param resolution The resolution of the output raster.
 * @param bounds The bounds of the raster.
 */
bool inRadius(double px, double py, int col, int row, double radius,
		double resolution, std::vector<double> &bounds){
	if(radius == 0.0) return true;
	// If a radius is given, extract the x and y of the current cell's centroid
	// and measure its distance (squared) from the point.
	double x = col * resolution + bounds[0] + resolution * 0.5;
	double y = row * resolution + bounds[1] + resolution * 0.5;
	// If the cell is outside the radius, ignore it.
	double r = sqrt(_sq(x - px) + _sq(y - py));
	return r <= radius;
}

void void_fill(Grid<float> &rast, double resolution) {
	int rlimit = _max(rast.rows(), rast.cols());
	double angle = 0.0;
	for(int r = 0; r < rast.rows(); ++r) {
		for(int c = 0; c < rast.cols(); ++c) {
			if(rast.get(c, r) != rast.nodata())
				continue;
			for(int rad = 1; rad < rlimit; ++rad) {
				// We'll average the values of pixels at this radius?
				float total = 0.0;
				int count = 0;
				for(double ang = 0; ang < 2.0 * PI; ang += std::atan(1.0 / rad)) {
					int r0 = (int) (r + std::sin(ang) * rad + 0.5);
					int c0 = (int) (c + std::cos(ang) * rad + 0.5);
					if(r0 < 0 || c0 < 0 || r0 >= rast.rows() || c0 >= rast.cols())
						continue;
					if(rast.get(c0, r0) != rast.nodata()) {
						total += rast.get(c0, r0);
						++count;
					}
				}
				if(count > 0) {
					rast.set(c, r, total / count);
					break;
				}
			}
		}
	}
}

void lasgrid(std::string &dstFile, std::vector<std::string> &files, std::set<int> &classes,
			int crs, int attribute, int type, double radius,
			double resolution, std::vector<double> &bounds, unsigned char angleLimit, bool fill) {

	if(files.size() == 0)
		_argerr("At least one input file is required.");
	_log(files.size() << " files.");

	if(dstFile.empty()) 
		_argerr("An output file is required.");
	_log("Creating " << dstFile);

	if(attribute == 0)
		_argerr("An attribute is required.");
	_log("Attribute: " << attribute);

	if(type == 0)
		_argerr("A valid type is required.");
	_log("Type: " << type);

	if(classes.size() == 0) {
		_log("WARNING: No classes given. Matching all classes.");
	} else {
		_log("Classes: " << classes.size());
	}

	if(angleLimit <= 0)
		_argerr("Angle limit must be greater than zero.");

	// Snap bounds
	if(bounds.size() == 4) 
		Util::snapBounds(bounds, resolution, 2);

	MemRaster<double> grid1;
	MemRaster<double> grid2;
	MemRaster<int> counts;
	MemRaster<std::vector<double>* > qGrid;

	las::ReaderFactory rf;
	std::vector<unsigned int> indices;

	bool hasBounds = bounds.size() > 0;
	if(!hasBounds) {
		bounds.push_back(DBL_MAX_POS);
		bounds.push_back(DBL_MAX_POS);
		bounds.push_back(DBL_MAX_NEG);
		bounds.push_back(DBL_MAX_NEG);
	}

	for(unsigned int i=0; i<files.size(); ++i) {

		_log("Checking file " << files[i]);

		std::ifstream in(files[i].c_str());
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();
		std::vector<double> bounds0 = { DBL_MAX_POS, DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_NEG };
		if(!Util::computeLasBounds(h, bounds0, 2))
			Util::computeLasBounds(r, bounds0, 2); // If the header bounds are bogus.
		in.close();
		if(hasBounds) {
			// If bounds are given, ignore blocks that do not intersect.
			if(!Util::intersects(bounds0, bounds, 2))
				continue;
		} else {
			// Otherwise expand to include each block.
			Util::expand(bounds, bounds0, 2);
		}
		indices.push_back(i);
	}

	Util::snapBounds(bounds, resolution, 2);
	Util::printBounds(bounds, 2);

	// Prepare grid
	int cols;
	int rows;
	Util::boundsToColsRows(bounds, resolution, &cols, &rows);

	// Compute the radius given the cell size, if it is not given.
	if(radius == -1.0)
		radius = sqrt(_sq(resolution / 2.0) * 2.0);

	_log("Raster size: " << cols << ", " << rows << "\nCell radius: " << radius);
	
	// For types other than count, we need a double grid to manage sums.
	if(type != TYPE_COUNT) {
		grid1.init(cols, rows);
		grid1.fill(-9999.0);
		// For variance and stddev, we need 2 double grids.
		if(type == TYPE_VARIANCE || type == TYPE_STDDEV || type == TYPE_PVARIANCE || type == TYPE_PSTDDEV) {
			grid2.init(cols, rows);
			grid2.fill(-9999.0);
		}
	}

	// For the median grid.
	if(type == TYPE_MEDIAN) {
		qGrid.init(cols, rows);
		qGrid.setDeallocator(*vector_dealloc);
		for(int i = 0; i < cols * rows; ++i)
			qGrid.set(i, new std::vector<double>());
	}

	// Create a grid for maintaining counts.
	counts.init(cols, rows);
	counts.fill(0);

	_log("Using " << indices.size() << " of " << files.size() << " files.");

	// Process files
	int current = 0; // Current file counter.
	for(unsigned int j = 0; j < indices.size(); ++j) {
		unsigned int i = indices[j];
		std::ifstream in(files[i].c_str());
		las::Reader reader = rf.CreateWithStream(in);
		las::Header header = reader.GetHeader();

		_log("File " << ++current << " of " << indices.size());

		while(reader.ReadNextPoint()) {
			las::Point pt = reader.GetPoint();
			
			if(_abs(pt.GetScanAngleRank()) > angleLimit)
				continue;

			double px = pt.GetX();
			double py = pt.GetY();
			// Check if in bounds, but only if clipping is desired.
			if(!Util::inBounds(px, py, bounds)) 
				continue;
			// If this point is not in the class list, skip it.
			if(!Util::inList(classes, pt.GetClassification().GetClass()))
				continue;
			// Get either the height or intensity value.
			double pz;
			if(attribute == ATT_INTENSITY) {
				pz = pt.GetIntensity();
			} else { // ATT_HEIGHT
				pz = pt.GetZ();
			}
			// Convert x and y, to col and row.
			int c = (int) ((px - bounds[0]) / resolution);
			int r = (int) ((py - bounds[1]) / resolution);
			// If the radius is > 0, compute the size of the window.
			int offset = radius > 0.0 ? (int) radius / resolution : 0;
			for(int cc = c - offset; cc < c + offset + 1; ++cc) {
				// Ignore out-of-bounds pixels.
				if(cc < 0 || cc >= cols) continue;
				for(int rr = r - offset; rr < r + offset + 1; ++rr) {
					// Ignore out-of-bounds pixels.
					if(rr < 0 || rr >= rows) continue;
					// If the coordinate is out of the cell's radius, continue.
					if(!inRadius(px, py, cc, rr, radius, resolution, bounds)) continue;
					// Compute the grid index. The rows are assigned from the bottom.
					int idx = (rows - rr - 1) * cols + cc;
					counts.set(idx, counts.get(idx) + 1);

					switch(type){
					case TYPE_MIN:
						if(counts[idx] == 1 || pz < grid1[idx])
							grid1.set(idx, pz);
						break;
					case TYPE_MAX:
						if(counts[idx] == 1 || pz > grid1[idx])
							grid1.set(idx, pz);
						break;
					case TYPE_VARIANCE:
					case TYPE_STDDEV:
					case TYPE_PVARIANCE:
					case TYPE_PSTDDEV:
						if(counts[idx] == 1) {
							grid2.set(idx, _sq(pz));
						} else {
							grid2.set(idx, grid2.get(idx) + _sq(pz));
						}
					case TYPE_MEAN:
						if(counts[idx] == 1) {
							grid1.set(idx, pz);
						} else {
							grid1.set(idx, grid1.get(idx) + pz);
						}
						break;
					case TYPE_QUANTILE:
					case TYPE_MEDIAN:
						qGrid[idx]->push_back(pz);
						break;
					}
				}
			}
		}
	}

	// Calculate cells or set nodata.
	switch(type) {
	case TYPE_MEAN:
		for(unsigned long i = 0; i < (unsigned long) cols * rows; ++i) {
			if(counts[i] > 0)
				grid1.set(i, grid1.get(i) / counts.get(i));
		}
		break;
	case TYPE_PVARIANCE:
		for(unsigned long i = 0; i < (unsigned long) cols * rows; ++i) {
			if(counts[i] > 1) {
				grid1.set(i, (grid2.get(i) - (_sq(grid1.get(i)) / counts.get(i))) / counts.get(i));
			} else {
				grid1.set(i, grid1.nodata());
			}
		}
		break;
	case TYPE_VARIANCE:
		for(unsigned long i = 0; i < (unsigned long) cols * rows; ++i) {
			if(counts[i] > 1) {
				grid1.set(i, (grid2.get(i) - (_sq(grid1.get(i)) / counts.get(i))) / (counts.get(i) - 1));
			} else {
				grid1.set(i, grid1.nodata());
			}
		}
		break;
	case TYPE_PSTDDEV:
		for(unsigned long i = 0; i < (unsigned long) cols * rows; ++i) {
			if(counts[i] > 1) {
				grid1.set(i, sqrt((grid2.get(i) - (_sq(grid1.get(i)) / counts.get(i))) / counts.get(i)));
			} else {
				grid1.set(i, grid1.nodata());
			}
		}
		break;
	case TYPE_STDDEV:
		for(unsigned long i = 0; i < (unsigned long) cols * rows; ++i) {
			if(counts[i] > 1) {
				grid1.set(i, sqrt((grid2.get(i) - (_sq(grid1.get(i)) / counts.get(i))) / (counts.get(i) - 1)));
			} else {
				grid1.set(i, grid1.nodata());
			}
		}
		break;
	case TYPE_DENSITY:
		{
			double r2 = _sq(resolution);
			for(unsigned long i = 0; i < (unsigned long) cols * rows; ++i) {
				if(counts[i] > 0) {
					grid1.set(i, counts[i] / r2);
				} else {
					grid1.set(i, 0.0);
				}
			}
		}
		break;
	case TYPE_MEDIAN:
		for(unsigned long i = 0; i < (unsigned long) cols * rows; ++i) {
			if(counts[i] > 0) {
				std::sort(qGrid[i]->begin(), qGrid[i]->end());
				int size = qGrid[i]->size();
				if(size % 2 == 0) {
					int idx = size / 2;
					grid1.set(i, ((*qGrid[i])[idx - 1] + (*qGrid[i])[idx]) / 2.0);
				} else {
					grid1.set(i, (*qGrid[i])[size / 2]);
				}
			}
		}
		break;
	}

	if(type == TYPE_COUNT) {
		Raster<int> rast(dstFile, bounds[0], bounds[1], bounds[2], bounds[3],
					resolution, -1, crs);
		rast.writeBlock(counts);
	} else {
		Raster<float> rast(dstFile, bounds[0], bounds[1], bounds[2], bounds[3],
					resolution, -9999.0, crs);
		// Cast the double grid to float for writing.
		MemRaster<float> tmp = (MemRaster<float>) grid1;
		tmp.nodata(rast.nodata());
		if(fill)
			void_fill(tmp, resolution);
		rast.writeBlock(tmp);
	}

}


int main(int argc, char **argv) {

	std::string dstFile;
	int crs = 0;
	int type = TYPE_MEAN;
	int att = ATT_HEIGHT;
	bool fill = false;
	double resolution = 2.0;
	double radius = -1.0;
	unsigned char angleLimit = 100;
	std::vector<double> bounds;
	std::set<int> classes;
	std::vector<std::string> files;

	for(int i = 1; i < argc; ++i) {
		std::string s(argv[i]);
		if(s == "-o") {
			dstFile = argv[++i];
		} else if(s == "-s") {
			crs = atoi(argv[++i]);
		} else if(s == "-f") {
			fill = true;
		} else if(s == "-t") {
			type = parseType(argv[++i]);
		} else if(s == "-r") {
			resolution = atof(argv[++i]);
		} else if(s == "-c") {
			Util::intSplit(classes, argv[++i]);
		} else if(s == "-a") {
			att = parseAtt(argv[++i]);
		} else if(s == "-d") {
			radius = atof(argv[++i]);
		} else if(s == "-v") {
			_loglevel(1);
		} else if(s == "--angle-limit") {
			angleLimit = (unsigned char) atoi(argv[++i]);
		} else if(s == "-b") {
			bounds.push_back(atof(argv[++i]));
			bounds.push_back(atof(argv[++i]));
			bounds.push_back(atof(argv[++i]));
			bounds.push_back(atof(argv[++i]));
		} else {
			files.push_back(argv[i]);
		}
	}

	try {
		lasgrid(dstFile, files, classes, crs, att, type, radius, resolution, bounds, angleLimit, fill);
	} catch(const std::exception &ex) {
		_log(ex.what());
		usage();
		return 1;
	} catch(const char *ex) {
		_log(ex);
		usage();
		return 1;
	}

	return 0;
}
