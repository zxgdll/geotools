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

#include "Util.hpp"
#include "Raster.hpp"

/**
 * Returns a value between 0 and 1 following the tan curve.
 * The range of step is expected to be 0 -> steps and is clamped.
 */
float _tanCurve(float step, float steps) {
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
void feather(Grid<float> &srcGrid, Grid<float> &dstGrid, int cols, int rows, float distance, float nodata, float resolution) {
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
					dstGrid.set(col, row, _tanCurve(step, steps)); // TODO: Configurable curves.
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
 * The distance argument determines the distance over which the edges are feathered. The
 * overviews argument, if true, forces the construction of new overviews. If this is not
 * used, the existing overviews may obscure the fact that the image has changed (when zoomed out.)
 */
 // TODO: Background must be larger than other files. Remedy this.
void mosaic(std::vector<std::string> &files, std::string &outfile, float distance) {
	
	if(distance <= 0.0)
		_argerr("Invalid distance: " << distance);

	if(outfile.size() == 0)
		_argerr("No output file given.");

	if(files.size() < 2)
		_argerr("Less than 2 files. Nothing to do.");

	_loglevel(1)

	// Copy the background file to the destination.
	_log("Copying background file.");
	Util::copyfile(files[0], outfile);

	// Open the destination file to modify it.
	Raster<float> output(outfile, 1, true);

	// Iterate over the files, adding each one to the background.
	for(unsigned int i = 1; i < files.size(); ++i) {

		_log("Processing file: " << files[i]);

		Raster<float> input(files[i]);
		if(output.resolutionX() != input.resolutionX() || output.resolutionY() != input.resolutionY()) 
			_argerr("Raster's (" << files[i] << ") resolution does not match the background (" << files[0] << ").");

		// Get the origin of the input w/r/t the output.
		int col = output.toCol(input.minx());
		int row = output.toRow(input.miny());
		int cols = output.toCol(input.maxx()) - col;
		int rows = output.toRow(input.maxy()) - row;

		float imNodata = input.nodata();
		float outNodata = output.nodata();

		// TODO: Configurable row height.
		int blkHeight = 500 > rows ? rows : 500; 					// The unbuffered height of a block
		int blkBuffer = (int) (distance / input.resolution()) + 1;  // The size of the buffer.

		_log("Loop rows: " << rows << ", " << blkHeight);

		// Move the window down by rowHeight and one rowOffset *not* two.
		#pragma omp parallel for
		for(int blkRow = 0; blkRow < (rows / blkHeight) + 1; blkRow += blkHeight) {

			// Initialize the grids.
			MemRaster<float> imGrid(cols, blkHeight + blkBuffer * 2);
			MemRaster<float> outGrid(cols, blkHeight + blkBuffer * 2);
			MemRaster<float> alphaGrid(cols, blkHeight + blkBuffer * 2);
			// Set to zero or nodata depending on the band.
			outGrid.fill(outNodata);
			imGrid.fill(imNodata);
			alphaGrid.fill(1.0);

			int blkRow0 = _max(0, blkRow - blkBuffer);
			int blkHeight0 = _min(rows - blkRow0 - 1, blkHeight + blkBuffer * 2);

			//std::cout << "Thread: " << omp_get_thread_num() << std::endl;
			_log("Indices: " << blkHeight << ", " << blkBuffer << ", " << blkRow << ", " << blkRow0 << ", " << blkHeight0);

			// Load the overlay.
			#pragma omp critical
			{
				imGrid.readBlock(0, blkRow0, cols, blkHeight0, imGrid);
			}

			_log("Feathering...");

			// Feather the overlay
			feather(imGrid, alphaGrid, cols, blkHeight0, distance, imNodata, input.resolution());

			// Read background data.
			#pragma omp critical 
			{
				output.readBlock(col, row + blkRow0, cols, blkHeight0, outGrid);
			}

			_log("Blending...");

			// Blend the overlay into the output.
			blend(imGrid, outGrid, alphaGrid, cols, blkHeight0, imNodata, outNodata);

			// Write back to the output.
			// We are extracting a slice out of the buffer, not writing the whole thing.
			#pragma omp critical 
			{
				output.writeBlock(col, row + blkRow0 + blkBuffer, cols, blkHeight0, outGrid);
			}

		}
	}
}

void usage() {
	std::cout << "Usage: mosaic [options] -o <output file> <file [file [file [...]]]>" << std::endl;
	std::cout << "    -o <file>     -- The output file." << std::endl;
	std::cout << "    -d <distance> -- The feather distance in map units (default 100.)" << std::endl;
	std::cout << "    <file [...]>  -- A list of files. The first is used as the background " << std::endl;
	std::cout << "                     and determines the size of the output. Subsequent " << std::endl;
	std::cout << "                     images are layered on top and feathered." << std::endl;
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

 		mosaic(files, outfile, distance);

 	} catch(const char *e) {
 		std::cerr << e << std::endl;
 		usage();
 		return 1;
 	}

 	return 0;
 }
