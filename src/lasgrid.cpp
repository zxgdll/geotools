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

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <liblas/liblas.hpp>

#include "Util.hpp"
#include "Raster.hpp"

#define TYPE_MIN 1
#define TYPE_MAX 2
#define TYPE_MEAN 3
#define TYPE_DENSITY 4
#define TYPE_VARIANCE 7
#define TYPE_STDDEV 8
#define TYPE_COUNT 9
#define TYPE_QUANTILE 10

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
	} else if(!strcmp("count", typeStr)) {
		return TYPE_COUNT;
	} else if(!strcmp("quantile", typeStr)) {
		return TYPE_QUANTILE;
	}
	return 0;
}

/**
 * Comparator for sorting floats.
 */
int _fcmp(const void * a, const void * b) {
	const float * aa = (const float *) a;
	const float * bb = (const float *) b;
	return (*aa > *bb) - (*bb > *aa);
}

void usage() {
	std::cerr << "Usage: lasgrid <options> <file [file [file]]>" << std::endl; 
	std::cerr << " -o <output file>" << std::endl;
	std::cerr << " -t <type>                   Output quantile, mean, max, min, variance, count, density, stddev (default mean)." << std::endl;
	std::cerr << " -r <resolution>             Resolution (default 2)." << std::endl;
	std::cerr << " -s <srid>                   The EPSG ID of the CRS." << std::endl;
	std::cerr << " -c <classes>                Comma-delimited (e.g. '2,0' (ground and unclassified))." << std::endl;
	std::cerr << " -a <attribute>              Use height, intensity (default height)." << std::endl;
	std::cerr << " -d <radius>                 Radius (not diameter); use zero for cell bounds." << std::endl;
	std::cerr << "                             For example, if the cell size is 2, the circumcircle's radius is sqrt(2) (1.414213562)." << std::endl;
	std::cerr << " -b <minx miny maxx maxy>    Extract points from the given box and create a raster of this size." << std::endl;
	std::cerr << " -q <num-quantiles,quantile> Gives the number of quantiles, and the index of the desired quantile. " << std::endl;
	std::cerr << "                             If there are n quantiles, there are n+1 possible indices, with 0 being " << std::endl;
	std::cerr << "                             the lower bound, and n+1 being the upper. For quartiles enter 4, for deciles 10, etc." << std::endl;
	std::cerr << " -v                          Verbose output." << std::endl;
}

void vector_dealloc(std::vector<float> *item) {
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
	if(radius == 0.0) return false;
	// If a radius is given, extract the x and y of the current cell's centroid
	// and measure its distance (squared) from the point.
	double x = col * resolution + bounds[0] + resolution * 0.5;
	double y = row * resolution + bounds[1] + resolution * 0.5;
	// If the cell is outside the radius, ignore it.
	double a = _sq(x - px) + _sq(y - py);
	double b = _sq(radius);
	return a <= b;
}

void lasgrid(std::string &dstFile, std::vector<std::string> &files, std::set<int> &classes,
			int crs, int quantile, int numQuantiles, int attribute, int type, float radius,
			float resolution, std::vector<double> &bounds, bool hasBounds, bool quiet) {

	// Snap bounds
	if(hasBounds) 
		Util::snapBounds(bounds, resolution, 2);

	Grid<float> grid1;
	Grid<float> grid2;
	Grid<int> counts;
	Grid<std::vector<float>* > qGrid;

	las::ReaderFactory rf;
	std::vector<unsigned int> indices;

	for(unsigned int i=0; i<files.size(); ++i) {
		std::ifstream in(files[i].c_str());
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();
		std::vector<double> bounds0 = { FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX };
		if(!Util::computeLasBounds(h, bounds0, 2))
			Util::computeLasBounds(r, bounds0, 2); // If the header bounds are bogus.
		in.close();
		// If bounds are given, ignore blocks that do not intersect.
		if(!hasBounds || (hasBounds && Util::intersects(bounds0, bounds, 2)))
			Util::expand(bounds, bounds0, 2);
		indices.push_back(i);
	}

	Util::snapBounds(bounds, resolution, 2);
	Util::printBounds(bounds, 2);

	// Prepare grid
	int cols = (int) ceil((bounds[2] - bounds[0]) / resolution);
	int rows = (int) ceil((bounds[3] - bounds[1]) / resolution);

	// Compute the radius given the cell size, if it is not given.
	if(radius == -1.0)
		radius = sqrt(_sq(resolution / 2.0) * 2.0);

	if(!quiet) {
		std::cerr << "Raster size: " << cols << ", " << rows << std::endl;
		std::cerr << "Cell radius: " << radius << std::endl;
	}
	
	// For types other than count, we need a float grid to manage sums.
	if(type != TYPE_COUNT) {
		grid1.init(cols, rows);
		grid1.fill(-9999.0);
		// For variance and stddev, we need 2 float grids.
		if(type == TYPE_VARIANCE || type == TYPE_STDDEV) {
			grid2.init(cols, rows);
			grid2.fill(-9999.0);
		}
	}

	// For the quantile grid.
	if(type == TYPE_QUANTILE) {
		qGrid.init(cols, rows);
		qGrid.setDeallocator(*vector_dealloc);
		for(int i = 0; i < cols * rows; ++i)
			qGrid[i] = new std::vector<float>();
	}

	// Create a grid for maintaining counts.
	counts.init(cols, rows);
	counts.fill(0);


	if(!quiet)
		std::cerr << "Using " << indices.size() << " of " << files.size() << " files." << std::endl;

	// Process files
	int current = 0; // Current file counter.
	for(unsigned int j = 0; j < indices.size(); ++j) {
		unsigned int i = indices[j];
		std::ifstream in(files[i].c_str());
		las::Reader reader = rf.CreateWithStream(in);
		las::Header header = reader.GetHeader();

		if(!quiet)
			std::cerr << "File " << ++current << " of " << indices.size() << std::endl;
	
		while(reader.ReadNextPoint()) {
			las::Point pt = reader.GetPoint();
			double px = pt.GetX();
			double py = pt.GetY();
			// Check if in bounds, but only if clipping is desired.
			if(hasBounds && !Util::inBounds(px, py, bounds)) 
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
			int offset = radius > 0.0 ? (int) ceil(radius / resolution) : 0;
			for(int cc = c - offset; cc < c + offset + 1; ++cc) {
				// Ignore out-of-bounds pixels.
				if(cc < 0 || cc >= cols) continue;
				for(int rr = r - offset; rr < r + offset + 1; ++rr) {
					// Ignore out-of-bounds pixels.
					if(rr < 0 || rr >= rows) continue;
					// If the coordinate is out of the cell's radius, continue.
					if(radius > 0.0 && !inRadius(px, py, cc, rr, radius, resolution, bounds)) continue;
					// Compute the grid index. The rows are assigned from the bottom.
					int idx = (rows - rr - 1) * cols + cc;
					switch(type){
					case TYPE_COUNT:
						counts[idx]++;
						break;
					case TYPE_MIN:
						if(counts[idx] == 0 || pz < grid1[idx])
							grid1[idx] = pz;
						counts[idx]++;
						break;
					case TYPE_MAX:
						if(counts[idx] == 0 || pz > grid1[idx])
							grid1[idx] = pz;
						counts[idx]++;
						break;
					case TYPE_MEAN:
					case TYPE_VARIANCE:
					case TYPE_STDDEV:
						if(counts[idx] == 0) {
							grid1[idx] = pz;
						} else {
							grid1[idx] += pz;
						}
						counts[idx]++;
						break;
					case TYPE_DENSITY:
						counts[idx]++;
						break;
					case TYPE_QUANTILE:
						qGrid[idx]->push_back(pz);
						counts[idx]++;
						break;
					}
				}
			}
		}
	}

	// Calculate cells or set nodata.
	switch(type) {
	case TYPE_MEAN:
	case TYPE_VARIANCE:
	case TYPE_STDDEV:
		for(unsigned long i = 0; i < (unsigned long) cols * rows; ++i) {
			if(counts[i] > 0)
				grid1[i] /= counts[i];
		}
		break;
	case TYPE_DENSITY:
		{
			double r2 = _sq(resolution);
			for(unsigned long i = 0; i < (unsigned long) cols * rows; ++i) {
				if(counts[i] > 0) {
					grid1[i] = counts[i] / r2;
				} else {
					grid1[i] = 0.0;
				}
			}
		}
		break;
	case TYPE_QUANTILE:
		for(unsigned long i = 0; i < (unsigned long) cols * rows; ++i) {
			if(counts[i] > numQuantiles) {
				// Compute the median.
				std::sort(qGrid[i]->begin(), qGrid[i]->end());
				if(quantile == 0) {
					// If index is zero, just return the min.
					grid1[i] = (*qGrid[i])[0];
				} else if(quantile == numQuantiles) {
					// If index == numQuantiles, return the max.
					grid1[i] = (*qGrid[i])[qGrid[i]->size() - 1];
				} else {
					float idx = (((float) qGrid[i]->size() - 1) / numQuantiles) * quantile;
					if(floor(idx) == idx) {
						grid1[i] = (*qGrid[i])[(int) idx];
					} else {
						grid1[i] = ((*qGrid[i])[(int) idx] + (*qGrid[i])[(int) idx + 1]) / 2.0;
					}
				}
			}
		}
		break;
	}

	// If we're doing the variance or stddev within cells, we have to run through the files again.
	if(type == TYPE_VARIANCE || type == TYPE_STDDEV) {
		current = 0;
		int c, r, cc, rr, idx;
		double px, py, pz;
		for(unsigned int j = 0; j < indices.size(); ++j) {
			unsigned int i = indices[j];
			std::ifstream in(files[i].c_str());
			las::Reader reader = rf.CreateWithStream(in);
			las::Header header = reader.GetHeader();

			if(!quiet) {
				std::cerr << "File " << ++current << " of " << indices.size() << std::endl;
				std::cerr << "Second run..." << std::endl;
			}	

			while(reader.ReadNextPoint()) {
				las::Point pt = reader.GetPoint();
				px = pt.GetX();
				py = pt.GetY();
				// Check if in bounds, but only if clipping is desired.
				if(hasBounds && !Util::inBounds(px, py, bounds)) 
					continue;
				// If point is not in class list, skip.
				if(!Util::inList(classes, pt.GetClassification().GetClass())) 
					continue;
				// Get value for height or intensity.
				if(attribute == ATT_INTENSITY) {
					pz = pt.GetIntensity();
				} else {
					pz = pt.GetZ();
				}

				// Translate point to row/col.
				c = (int) ((px - bounds[0]) / resolution);
				r = (int) ((py - bounds[1]) / resolution);
				// If the radius is > 0, compute the cell offset to accomodate it.
				int offset = radius > 0.0 ? (int) ceil(radius / resolution) : 0;
				for(cc = c - offset; cc < c + offset + 1; ++cc) {
					// Skip if out of bounds.
					if(cc < 0 || cc >= cols) 
						continue;
					for(rr = r - offset; rr < r + offset + 1; ++rr) {
						// Skip if out of bounds.
						if(rr < 0 || rr >= rows) 
							continue;
						// If the coordinate is out of the cell's radius, continue.
						if(radius > 0.0 && !inRadius(px, py, cc, rr, radius, resolution, bounds))
							continue;
						// The rows are assigned from the bottom.
						idx = (rows - rr - 1) * cols + cc;
						switch(type){
						case TYPE_VARIANCE:
						case TYPE_STDDEV:
							grid2[idx] += _sq(grid1[idx] - pz);
							break;
						}
					}
				}
			}
		}

		// Calculate cells or set nodata.
		switch(type) {
		case TYPE_VARIANCE:
			for(int i=0;i<cols*rows;++i) {
				if(counts[i] > 0)
					grid1[i] = grid2[i] / counts[i];
			}
			break;
		case TYPE_STDDEV:
			for(int i=0;i<cols*rows;++i) {
				if(counts[i] > 0)
					grid1[i] = sqrt(grid2[i] / counts[i]);
			}
			break;
		}
	}

	Raster<float> rast(dstFile, bounds[0], bounds[1], bounds[2], bounds[3],
				resolution, -9999.0, crs);
	rast.writeBlock(grid1);

}

int main(int argc, char **argv) {

	std::string dstFile;
	bool hasBounds = false;
	bool quiet = true;
	int crs = 0;
	int type = TYPE_MEAN;
	int att = ATT_HEIGHT;
	int numQuantiles = 0;
	int quantile = 0;
	float resolution = 2.0;
	float radius = -1.0;
	std::vector<double> bounds = {FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX};
	std::set<int> classes;
	std::vector<int> quantiles;
	std::vector<std::string> files;

	for(int i = 1; i < argc; ++i) {
		std::string s(argv[i]);
		if(s == "-o") {
			dstFile = argv[++i];
		} else if(s == "-s") {
			crs = atoi(argv[++i]);
		} else if(s == "-t") {
			type = parseType(argv[++i]);
		} else if(s == "-r") {
			resolution = atof(argv[++i]);
		} else if(s == "-c") {
			Util::intSplit(classes, argv[++i]);
		} else if(s == "-q") {
			Util::intSplit(quantiles, argv[++i]);
		} else if(s == "-a") {
			att = parseAtt(argv[++i]);
		} else if(s == "-d") {
			radius = atof(argv[++i]);
		} else if(s == "-v") {
			quiet = false;
		} else if(s == "-b") {
			bounds[0] = atof(argv[++i]);
			bounds[1] = atof(argv[++i]);
			bounds[2] = atof(argv[++i]);
			bounds[3] = atof(argv[++i]);
			hasBounds = true;
		} else {
			files.push_back(argv[i]);
		}
	}

	if(files.size() == 0) {
		std::cerr << "At least one input file is required." << std::endl;
		usage();
		return 1;
	} else {
		std::cerr << files.size() << " files." << std::endl;
	}

	if(dstFile.empty()) {
		std::cerr << "An output file is required." << std::endl;
		usage();
		return 1;
	} else {
		if(!quiet)
			std::cerr << "Creating " << dstFile << std::endl;
	}

	if(att == 0) {
		std::cerr << "An attribute is required." << std::endl;
		usage();
		return 1;
	} else {
		if(!quiet)
			std::cerr << "Attribute: " << att << std::endl;
	}

	if(type == 0) {
		std::cerr << "A valid type is required." << std::endl;
		usage();
		return 1;
	} else {
		if(!quiet)
			std::cerr << "Type: " << type << std::endl;
	}

	if(type == TYPE_QUANTILE && quantiles.size() < 2) {
		std::cerr << "If the type is quantiles, there must be a -q argument with 2 values." << std::endl;
		usage();
		return 1;
	} else if(type == TYPE_QUANTILE) {
		quantile = quantiles[1];
		numQuantiles = quantiles[0];
		if(quantile == 0) {
			if(!quiet)
				std::cerr << "Using min because q = 0." << std::endl;
			type = TYPE_MIN;
		} else if(quantile == numQuantiles) {
			if(!quiet)
				std::cerr << "Using max because q = numQuantiles." << std::endl;
			type = TYPE_MAX;
		} else {
			if(numQuantiles < 2) {
				std::cerr << "The number of quantiles must be >=2 and the quantile must be between 0 and n, inclusive." << std::endl;
				usage();
				return 1;
			}
			if(!quiet)
				std::cerr << "Quantiles: " << quantile << ", Q: " << numQuantiles << std::endl;
		}
	}

	if(classes.size() == 0) {
		std::cerr << "WARNING: No classes given. Matching all classes." << std::endl;
	} else {
		if(!quiet)
			std::cerr << "Classes: " << classes.size() << std::endl;
	}

	try {
		lasgrid(dstFile, files, classes, crs, quantile, numQuantiles, att,
				type, radius, resolution, bounds, hasBounds, quiet);
	} catch(const char *ex) {
		std::cerr << ex << std::endl;
		return 1;
	}
}
