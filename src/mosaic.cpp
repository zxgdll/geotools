/**
 * "Feathers" the edges of a data region (i.e. not nodata) to the specified
 * distance in map units, using the specified curve.
 */

//#include <float.h>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#define PI 3.14159265358979323846
#define E 0.5772156649
#define SLOPE 56.7

float sinTable[8];
float cosTable[8];

/**
 * Initialize the sin and cos lookup tables. This way we don't
 * have to compute them in each run, nor use angles.
 */
void initSinCos() {
	int i = 0;
	for(float d = 0.0; d <= PI * 2.0; d += PI / 4.0, ++i) {
		sinTable[i] = sin(d);
		cosTable[i] = cos(d);
	}
}

/**
 * Returns the minimum distance of the given cell from a null cell as 
 * distinguished by nodata. The limit distance, beyond which the
 * algorithm will not search, is in map units. The algorithm searches
 * in a circular pattern, checking for nulls at each PI/4 interval.
 * I if it finds one, it returns the distance. If not, it continues searching
 * further out.
 */
float distanceFromNull(float *srcGrid, int cols, int rows, int col, int row, float nodata, float resolution, float limit) {
	if(srcGrid[row * cols + col] == nodata)
		return 0.0;
	int r, c;
	for(float d = 0.0; d <= limit; d += resolution) {
		for(int a = 0; a < 8; ++a) {
			c = col + (int) round(cosTable[a] * d);
			r = row + (int) round(sinTable[a] * d);
			if(c >= 0 && c < cols && r >= 0 && r < rows
				&& srcGrid[r * cols + c] == nodata) {
				return d;
			}
		}
	}
	return limit;
}

/**
 * Compute a logistic curve based on the distance (variable) and the
 * limit (static).
 */
inline float logCurve(float distance, float limit) {
	return 1.0 - 1.0 / (1.0 + pow(E, -(1.0 / limit * SLOPE) * (distance - limit / 2.0)));
}

/**
 * Feathers the edges of data regions in a raster grid by setting the alpha
 * value for a pixel in proportion to its distance from the nearest null.
 */
void feather(float *srcGrid, float *dstGrid, int cols, int rows, float distance, float nodata, float resolution) {
	float dfn = 0.0;
	for(int row = 0; row < rows; ++row) {
		for(int col = 0; col < cols; ++col) {
			dfn = distanceFromNull(srcGrid, cols, rows, col, row, nodata, resolution, distance);
			dstGrid[row * cols + col] =  logCurve(dfn, distance); //dfn / distance;
		}
	}
}

/**
 * Blends two rasters together using the alpha grid for blending.
 */
void blend(float *imgGrid, float *bgGrid, float *alpha, int cols, int rows, float imNodata, float bgNodata) {
	// TODO: Maybe im component needs to be scaled down by the number of images it's blending with?
	for(int i = 0; i < cols * rows; ++i) {
		if(bgGrid[i] == bgNodata || imgGrid[i] == imNodata)
			continue;
		bgGrid[i] = bgGrid[i] * (1.0 - alpha[i]) + imgGrid[i] * alpha[i];
	}
}


int build(std::vector<std::string> &files, std::string &outfile, float distance) {
	if(files.size() < 2) {
		std::cerr << "Less than 2 files. Nothing to do." << std::endl;
		return 1;
	}

	GDALAllRegister();

	// Copy the background file to the destination.
	std::ifstream src(files[0].c_str(), std::ios::binary);
	std::ofstream dst(outfile.c_str(), std::ios::binary);
	dst << src.rdbuf();

	GDALDataset *outDS = NULL, *imDS = NULL;
	float *imGrid = NULL, *outGrid = NULL, *alphaGrid = NULL, *tmpGrid = NULL;
	double outTrans[6], imTrans[6];
	float imNodata, outNodata;
	int celBuf;
	int cols, rows, imCols, imRows;
	int col, row;

	try {

		// Open the destination file to modify it.
		outDS = (GDALDataset *) GDALOpen(outfile.c_str(), GA_Update);
		if(outDS == NULL)
			throw "Failed to open destination image.";

		// Load the transform parameters.
		outDS->GetGeoTransform(outTrans);

		// Compute the cell buffer required to accomodate the distance argument.
		celBuf = (int) (distance / outTrans[1]) + 1;

		// Iterate over the files, adding each one to the background.
		for(unsigned int i = 1; i < files.size(); ++i) {

			std::cout << "Processing file: " << files[i] << std::endl;

			imDS = (GDALDataset *) GDALOpen(files[i].c_str(), GA_ReadOnly);
			if(!imDS)
				throw "Failed to open mosaic image.";

			imDS->GetGeoTransform(imTrans);
			if(outTrans[1] != imTrans[1] || outTrans[5] != imTrans[5]) {
				GDALClose(imDS);
				throw "Raster's resolution does not match the background.";
			}

			cols = imDS->GetRasterXSize();
			rows = imDS->GetRasterYSize();
			col = (int) ((imTrans[0] - outTrans[0]) / outTrans[1]);
			row = (int) ((imTrans[3] - outTrans[3]) / outTrans[5]);
			imNodata = imDS->GetRasterBand(1)->GetNoDataValue();
			imCols = cols + celBuf * 2;
			imRows = rows + celBuf * 2;
			outNodata = outDS->GetRasterBand(1)->GetNoDataValue();
			
			// Initialize the grids.
			imGrid = (float *) malloc(imCols * imRows * sizeof(float));
			alphaGrid = (float *) calloc(imCols * imRows, sizeof(float));

			// Fill the im grid with nodata.
			std::fill_n(imGrid, imCols * imRows, imNodata);

			// Load the overlay into a tmp grid, then copy it into
			// the im grid with the correct buffer.
			tmpGrid = (float *) malloc(cols * rows * sizeof(float));
			imDS->RasterIO(GF_Read, 0, 0, cols, rows, tmpGrid, cols, rows, GDT_Float32, 1, NULL, 0, 0, 0, NULL);
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c) {
					imGrid[(r + celBuf) * imCols + (c + celBuf)] = tmpGrid[r * cols + c];
				}
			}
			free(tmpGrid);

			// Initialize an out grid of the correct size for the feathered overlay.			
			outGrid = (float *) malloc(imCols * imRows * sizeof(float));
			outDS->RasterIO(GF_Read, col - celBuf, row - celBuf, imCols, imRows, outGrid, imCols, imRows, GDT_Float32, 1, NULL, 0, 0, 0, NULL);

			std::cout << "Feathering" << std::endl;
			// Feather the overlay
			feather(imGrid, alphaGrid, imCols, imRows, distance, imNodata, outTrans[1]);

			std::cout << "Blending" << std::endl;
			// Blend the overlay into the output.
			blend(imGrid, outGrid, alphaGrid, imCols, imRows, imNodata, outNodata);

			// Write back to the output.
			outDS->RasterIO(GF_Write, col - celBuf, row - celBuf, imCols, imRows, outGrid, imCols, imRows, GDT_Float32, 1, NULL, 0, 0, 0, NULL);

			free(imGrid);
			free(outGrid);
			free(alphaGrid);

			GDALClose(imDS);
		}
	} catch(const char *e) {

		free(imGrid);
		free(outGrid);
		free(alphaGrid);
		free(tmpGrid);

		GDALClose(outDS);
		GDALClose(imDS);
		throw e;
	}

	GDALClose(outDS);

	return 0;
}

int main(int argc, char **argv) {

 	try {

	 	float distance = 100.0;
	 	std::vector<std::string> files;
	 	std::string outfile;

	 	for(int i = 1; i < argc; ++i) {
	 		std::string arg(argv[i]);
	 		if(arg == "-d") {
	 			distance = atof(argv[++i]);
	 		} else if(arg == "-o") {
	 			outfile = argv[++i];
	 		} else {
	 			files.push_back(argv[i]);
	 		}
	 	}

	 	initSinCos();

 		return build(files, outfile, distance);

 	} catch(const char *e) {
 		std::cerr << e << std::endl;
 		return 1;
 	}

 }
