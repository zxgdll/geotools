/**
 * "Feathers" the edges of a data region (i.e. not nodata) to the specified
 * distance in map units, using the specified curve.
 */
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <omp.h>

#include "geotools.h"
#include "Util.hpp"
#include "Raster.hpp"

using namespace raster;

/**
 * Returns a value between 0 and 1 following the tan curve.
 * The range of step is expected to be 0 -> steps and is clamped.
 */
float tanCurve(float step, float steps) {
	step = _min(steps, _max(0.0, step));
	return tanh(((step - steps / 2.0) / (steps / 2.0)) * PI) * 0.5 + 0.5;
}

/**
 * Returns true if the given pixel is next to a nodata pixel (but not if
 * it is one), or the edge of the grid.
 */
bool isEdgePixel(Grid<char> &fillGrid, int col, int row, int cols, int rows) {
	if(fillGrid.get(col, row) == 0)
		return false;
	for(int r = row - 1; r < row + 2; ++r) {
		for(int c = col - 1; c < col + 2; ++c) {
			if(c <= 0 || r <= 0 || c >= cols-1 || r >= rows-1 || fillGrid.get(c, r) == 0)
				return true;
		}
	}
	return false;
}

/**
 * Feathers the edges of data regions in a raster grid by setting the alpha
 * value for a pixel in proportion to its distance from the nearest null or edge.
 */
void feather(Grid<float> &srcGrid, Grid<float> &dstGrid, int cols, int rows, float distance, float nodata, double resolution) {
	// Fill grid is used to keep track of where the edges are as they're "snowed in"
	// Starts out as a mask of non-nodata pixels from the source.
	MemRaster<char> fillGrid(cols, rows);
	for(unsigned long i = 0; i < (unsigned long) rows * cols; ++i)
		fillGrid.set(i, srcGrid[i] == nodata ? 0 : 1);
	// The number of steps is just the number of pixels needed to 
	// cover the distance of the fade.
	float step = 0.0;
	float steps = _max(1.0, distance / resolution);
	bool found = false;
	// "Snow in" the alpha mask. The loop will exit when no more edge pixels can be found.
	do {
		found = false;
		for(int row = 0; row < rows; ++row) {
			for(int col = 0; col < cols; ++col) {
				if(isEdgePixel(fillGrid, col, row, cols, rows)) {
					fillGrid.set(col, row, 2); // Set edge to dirty.
					dstGrid.set(col, row, tanCurve(step, steps)); // TODO: Configurable curves.
					found = true;
				}
			}
		}
		// Reset dirty edges to 0
		for(int i = 0; i < rows * cols; ++i) 
			if(fillGrid[i] == 2) fillGrid.set(i, 0);
		++step;
	} while(found);
}

/**
 * Blends two rasters together using the alpha grid for blending.
 */
void blend(Grid<float> &imgGrid, Grid<float> &bgGrid, Grid<float> &alpha, int cols, int rows, float imNodata, float bgNodata) {
	for(unsigned long i = 0; i < (unsigned long) cols * rows; ++i) {
		if(!(bgGrid[i] == bgNodata || imgGrid[i] == imNodata))
			bgGrid.set(i, bgGrid[i] * (1.0 - alpha[i]) + imgGrid[i] * alpha[i]);
	}
}

/**
 * Mosaic the given files together using the first as the base. The base file will serve as a clipping
 * mask for the others. The spatial reference system and resolution of all layers must match.
 * The distance argument determines the distance over which the edges are feathered. 
 */
 // TODO: Background must be larger than other files. Remedy this.
void mosaic(std::vector<std::string> &files, std::string &outfile, float distance, int rowHeight) {
	
	if(distance <= 0.0)
		_argerr("Invalid distance: " << distance);
	if(outfile.size() == 0)
		_argerr("No output file given.");
	if(files.size() < 2)
		_argerr("Less than 2 files. Nothing to do.");
	if(rowHeight < 1)
		_argerr("Less than 1 row given.");

	// Copy the background file to the destination.
	_trace("Copying background file.");
	Util::copyfile(files[0], outfile);

	// Open the destination file to modify it.
	Raster<float> output(outfile, 1, true);

	// Iterate over the files, adding each one to the background.
	for(unsigned int i = 1; i < files.size(); ++i) {

		_trace("Processing file: " << files[i]);

		Raster<float> input(files[i]);
		if(output.resolutionX() != input.resolutionX() || output.resolutionY() != input.resolutionY()) 
			_argerr("Raster's (" << files[i] << ") resolution does not match the background (" << files[0] << ").");

		// Get the origin of the input w/r/t the output.
		int col = output.toCol(input.lx());
		int row = output.toRow(input.ty());
		int cols = output.toCol(input.rx()) - col;
		int rows = output.toRow(input.by()) - row;

		float imNodata = input.nodata();
		float outNodata = output.nodata();

		int blkHeight = _min(rows, rowHeight);				        // The unbuffered height of a block
		int blkBuffer = (int) (distance / input.resolutionX()) + 1; // The size of the buffer.

		// Initialize the grids.
		MemRaster<float> imGrid(cols, blkHeight + blkBuffer * 2);
		MemRaster<float> outGrid(cols, blkHeight + blkBuffer * 2);
		MemRaster<float> alphaGrid(cols, blkHeight + blkBuffer * 2);

		#pragma omp parallel for
		for(int blkRow = blkBuffer; blkRow < input.rows(); blkRow += blkHeight) {

			// Set to zero or nodata depending on the band.
			outGrid.fill(outNodata);
			imGrid.fill(imNodata);
			alphaGrid.fill(1.0);

			int blkRow0 = blkRow - blkBuffer;
			int blkHeight0 = _min(input.rows() - blkRow0 - 1, blkHeight + blkBuffer * 2);

			_trace("Reading source block...")
			// Load the overlay.
			#pragma omp critical
			{
				input.readBlock(0, blkRow0, cols, blkHeight0, imGrid);
			}

			_trace("Feathering...")
			// Feather the overlay
			feather(imGrid, alphaGrid, cols, blkHeight0, distance, imNodata, input.resolutionX());

			// Read background data.
			#pragma omp critical 
			{
				output.readBlock(col, row + blkRow0, cols, blkHeight0, outGrid);
			}

			_trace("Blending...")
			// Blend the overlay into the output.
			blend(imGrid, outGrid, alphaGrid, cols, blkHeight0, imNodata, outNodata);

			_trace("Writing to output...")
			// Write back to the output.
			// We are extracting a slice out of the buffer, not writing the whole thing.
			#pragma omp critical 
			{
				std::unique_ptr<Grid<float> > slice = outGrid.slice(0, blkBuffer, cols, blkHeight);
				output.writeBlock(col, row + blkRow, cols, blkHeight, *slice);
			}

		}
	}
}

void usage() {
	std::cerr << "Usage: mosaic [options] -o <output file> <file [file [file [...]]]>\n"
		<< "    -o <file>     -- The output file.\n"
		<< "    -d <distance> -- The feather distance in map units (default 100.)\n"
		<< "    -b <rows>     -- The height of working block, in pixels. The distance buffer is added to this.\n"
		<< "                     default is 500. A smaller or larger value can be used to conserve or exploit\n"
		<< "                     memory.\n"
		<< "    <file [...]>  -- A list of files. The first is used as the background\n"
		<< "                     and determines the size of the output. Subsequent \n"
		<< "                     images are layered on top and feathered.\n"
		<< "    -v               Verbose mode.\n";
}

int main(int argc, char **argv) {

 	try {

	 	float distance = 100.0;
	 	int rows = 500;
	 	std::vector<std::string> files;
	 	std::string outfile;
	 	bool verbose = false;

	 	for(int i = 1; i < argc; ++i) {
	 		std::string arg(argv[i]);
	 		if(arg == "-d") {
	 			distance = atof(argv[++i]);
	 		} else if(arg == "-b") {
	 			rows = atoi(argv[++i]);
	 		} else if(arg == "-o") {
	 			outfile = argv[++i];
	 		} else if(arg == "-v") {
	 			verbose = true;
	 		} else {
	 			files.push_back(argv[i]);
	 		}
	 	}

	 	_loglevel(verbose ? LOG_TRACE : LOG_ERROR);

 		mosaic(files, outfile, distance, rows);

 	} catch(const std::exception &e) {
 		std::cerr << e.what() << std::endl;
 		usage();
 		return 1;
 	}

 	return 0;
 }
