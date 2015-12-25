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

//#include <omp.h>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <liblas/liblas.hpp>

#include "Util.hpp"

#define TYPE_MIN 1
#define TYPE_MAX 2
#define TYPE_MEAN 3
#define TYPE_DENSITY 4
#define TYPE_VARIANCE 7
#define TYPE_STDDEV 8
#define TYPE_COUNT 9

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
	}
	return 0;
}

/**
 * Square a float.
 */
inline float _sq(float a) {
	return a * a;
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
	std::cerr << " -t <type>                   mean, max, min, variance, count, density, stddev (default mean)" << std::endl;
	std::cerr << " -r <resolution>             (default 2)" << std::endl;
	std::cerr << " -s <srid>                   The EPSG ID of the CRS." << std::endl;
	std::cerr << " -c <classes>                comma-delimited (e.g. '2,0' (ground and unclassified))" << std::endl;
	std::cerr << " -a <attribute>              height, intensity (default height)" << std::endl;	
	std::cerr << " -d <radius>                 use zero for cell bounds" << std::endl;	
	std::cerr << " -b <minx miny maxx maxy>    extract points from the given box and create a raster of this size" << std::endl;	
}

int main(int argc, char **argv) {

	char *dstFile = NULL;
	int crs = 0;
	int type = TYPE_MEAN;
	int att = ATT_HEIGHT;
	double resolution = 2.0;
	double radius = -1.0;
	std::set<int> classes;
	std::vector<std::string> files;
	double bounds[] = {FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX};
	bool hasBounds = false;
	bool quiet = false;
	
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
		} else if(s == "-a") {
			att = parseAtt(argv[++i]);
		} else if(s == "-d") {
			radius = atof(argv[++i]);
		} else if(s == "-q") {
			quiet = true;
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

	if(dstFile == NULL) {
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

	if(classes.size() == 0) {
		std::cerr << "WARNING: No classes given. Matching all classes." << std::endl;
	} else {
		if(!quiet)
			std::cerr << "Classes: " << classes.size() << std::endl;
	}

	// Snap bounds
	if(hasBounds) 
		Util::snapBounds(bounds, resolution, 2);

	int cols, rows;
	float *grid = NULL, *grid2 = NULL;
	int *counts = NULL;

	las::ReaderFactory rf;
	std::vector<unsigned int> indices;

	for(unsigned int i=0; i<files.size(); ++i) {
		std::ifstream in(files[i].c_str());
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();
		double bounds0[] = { FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX };
		Util::computeLasBounds(h, bounds0, 2);
		in.close();
		// If bounds are given, ignore blocks that do not intersect.
		if(hasBounds) {
			if(!Util::intersects(bounds0, bounds, 2))
				continue;
		} else {
			Util::expand(bounds, bounds0, 2);
		}
		indices.push_back(i);
	}

	Util::snapBounds(bounds, resolution, 2);
	Util::printBounds(bounds, 2);

	// Prepare grid
	cols = (int) ceil((bounds[2] - bounds[0]) / resolution);
	rows = (int) ceil((bounds[3] - bounds[1]) / resolution);
	std::cerr << "Raster size: " << cols << ", " << rows << std::endl;

	// Compute the radius given the cell size, if it is not given.
	if(radius == -1.0)
		radius = sqrt(_sq(resolution / 2.0) * 2);
	if(!quiet)
		std::cerr << "Cell radius: " << radius << std::endl;
	
	// For types other than count, we need a double grid to manage sums.
	if(type != TYPE_COUNT) {
		grid = (float *) calloc(sizeof(float), cols * rows);
		if(grid == NULL) {
			std::cerr << "Failed to initialize grid." << std::endl;
			return 1;
		}

		// For variance and stddev, we need 2 double grids.
		if(type == TYPE_VARIANCE || type == TYPE_STDDEV) {
			grid2 = (float *) calloc(sizeof(float), cols * rows);
			if(grid2 == NULL) {
				std::cerr << "Failed to initialize secondary grid." << std::endl;
				return 1;
			}
		}
	}

	// Create a grid for maintaining counts.
	counts = (int *) calloc(sizeof(int), cols * rows);
	if(counts == NULL) {
		std::cerr << "Failed to initialize counts." << std::endl;
		return 1;
	}

	int current = 0;

	if(!quiet)
		std::cerr << "Using " << indices.size() << " of " << files.size() << " files." << std::endl;
	// Process files
	int c, r, cc, rr, idx;
	double x, y;
	double px, py, pz;
	for(unsigned int j = 0; j < indices.size(); ++j) {
		unsigned int i = indices[j];
		std::ifstream in(files[i].c_str());
		las::Reader reader = rf.CreateWithStream(in);
		las::Header header = reader.GetHeader();

		if(!quiet)
			std::cerr << "File " << ++current << " of " << indices.size() << std::endl;
	
		while(reader.ReadNextPoint()) {
			las::Point pt = reader.GetPoint();
			px = pt.GetX();
			py = pt.GetY();
			// Check if in bounds, but only if clipping is desired.
			if(hasBounds && !Util::inBounds(px, py, bounds)) 
				continue;
			// If this point is not in the class list, skip it.
			if(!Util::inList(classes, pt.GetClassification().GetClass()))
				continue;
			// Get either the height or intensity value.
			if(att == ATT_INTENSITY) {
				pz = pt.GetIntensity();
			} else { // ATT_HEIGHT
				pz = pt.GetZ();
			}
			// Convert x and y, to col and row.
			c = (int) ((px - bounds[0]) / resolution);
			r = (int) ((py - bounds[1]) / resolution);

			// If the radius is > 0, compute the size of the window.
			int offset = radius > 0.0 ? (int) ceil(radius / resolution) : 0;
			for(cc = c - offset; cc < c + offset + 1; ++cc) {
				// Ignore out-of-bounds pixels.
				if(cc < 0 || cc >= cols) 
					continue;
				for(rr = r - offset; rr < r + offset + 1; ++rr) {
					// Ignore out-of-bounds pixels.
					if(rr < 0 || rr >= rows) 
						continue;
					// If a radius is given, extract the x and y of the current cell's centroid
					// and measure its distance (squared) from the point.
					if(radius > 0.0) {
						// TODO: In UTM N, the resolution of the vertical is negative
						// because the largest coordinate is at the top. We pretend that
						// all projections are UTMN.
						x = cc * resolution + bounds[0] + resolution * 0.5; 
						y = rr * resolution + bounds[1] + resolution * 0.5;
						// If the cell is outside the radius, ignore it.
						if((_sq(x - px) + _sq(y - py)) > (radius * radius))
							continue;
					}
					// Compute the grid index. The rows are assigned from the bottom.
					idx = (rows - rr - 1) * cols + cc;
					switch(type){
					case TYPE_COUNT:
						counts[idx]++;
						break;
					case TYPE_MIN:
						if(counts[idx] == 0 || pz < grid[idx]) 
							grid[idx] = pz;
						counts[idx]++;
						break;
					case TYPE_MAX:
						if(counts[idx] == 0 || pz > grid[idx])
							grid[idx] = pz;
						counts[idx]++;
						break;
					case TYPE_MEAN:
					case TYPE_VARIANCE:
					case TYPE_STDDEV:
						grid[idx] += pz;
						counts[idx]++;
						break;
					case TYPE_DENSITY:
						grid[idx]++;
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
		for(int i=0;i<cols*rows;++i) {
			if(counts[i] > 0) {
				grid[i] /= counts[i];
			} else {
				grid[i] = -9999.0;
			}
		}
		break;
	case TYPE_DENSITY:
		{
		double r2 = _sq(resolution);
		for(int i=0;i<cols*rows;++i)
			if(counts[i] > 0) {
				grid[i] /= r2;
			} else {
				grid[i] = 0.0;
			}
		}	
		break;
	case TYPE_COUNT:
		// do nothing
		break;
	default:
		for(int i=0;i<cols*rows;++i) {
			if(counts[i] == 0)
				grid[i] = -9999.0;
		}
		break;
	}

	// If we're doing the variance or stddev within cells, we have to run through the files again.
	if(type == TYPE_VARIANCE || type == TYPE_STDDEV) {
		current = 0;
		int c, r, cc, rr, idx;
		double x, y, px, py, pz;
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
				if(att == ATT_INTENSITY) {
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
						// If a radius is given compute the position of the cell centroid.
						if(radius > 0.0) {
							x = cc * resolution + bounds[0] + resolution * 0.5; 
							y = rr * resolution + bounds[1] + resolution * 0.5;
							// If the point is outside of the radius, skip it.
							if((_sq(x - px) + _sq(y - py)) > radius * radius) 
								continue;
						}
						// The rows are assigned from the bottom.
						idx = (rows - rr - 1) * cols + cc;
						switch(type){
						case TYPE_VARIANCE:
						case TYPE_STDDEV:
							grid2[idx] += _sq(grid[idx] - pz);
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
					grid[i] = grid2[i] / counts[i];
			}
			break;
		case TYPE_STDDEV:
			for(int i=0;i<cols*rows;++i) {
				if(counts[i] > 0)
					grid[i] = sqrt(grid2[i] / counts[i]);
			}
			break;
		}

		free(grid2);
	}

	GDALAllRegister();

	// Save raster
	const char *pszFormat = "GTiff";
	char *wkt = NULL;
	char **papszOptions = NULL;
	double geoTransform[] = {bounds[0], resolution, 0.0, bounds[3], 0.0, -resolution};
	GDALDataType gType = GDT_Float32;
	GDALDriver *poDriver;
	GDALDataset *poDstDS;       
	OGRSpatialReference ref;
	GDALRasterBand *band;
	void *data;

	if(type == TYPE_COUNT) {
		gType = GDT_Int32;
		data = (void *) counts;
	} else {
		data = (void *) grid;
	}

   	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	poDstDS = poDriver->Create( dstFile, cols, rows, 1, gType, papszOptions );
	poDstDS->SetGeoTransform(geoTransform);
    	
	if(crs > 0) {
		ref.importFromEPSG(crs);
		ref.exportToWkt(&wkt);
   		poDstDS->SetProjection(wkt);
	}

	band = poDstDS->GetRasterBand(1);
	band->SetNoDataValue(-9999.0);
	CPLErr err = band->RasterIO(GF_Write, 0, 0, cols, rows, data, cols, rows, gType, 0, 0);
	if(0 != err)
		std::cerr << "Failed to write raster.";
	GDALClose(poDstDS);
	CPLFree(wkt);
	return err; // TODO: Different return method.

}
