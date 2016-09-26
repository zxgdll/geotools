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
#include "util.hpp"
#include "raster.hpp"

#include "raster.cpp"

using namespace geotools::raster;

namespace geotools {

	namespace raster {

		namespace util {

			/**
			 * Returns a value between 0 and 1 following the tan curve.
			 * The range of step is expected to be 0 -> steps and is clamped.
			 */
			float tanCurve(float step, float steps) {
				step = g_min(steps, g_max(0.0, step));
				return tanh(((step - steps / 2.0) / (steps / 2.0)) * G_PI) * 0.5 + 0.5;
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
				g_debug("feather: nodata: " << nodata << "; res: " << resolution << "; distance: " << distance << "; cols: " << cols << "; rows: " << rows);
				MemRaster<char> fillGrid(cols, rows);
				for(size_t i = 0; i < (size_t) rows * cols; ++i)
					fillGrid.set(i, srcGrid[i] == nodata ? 0 : 1);

				// The number of steps is just the number of pixels needed to 
				// cover the distance of the fade.
				float step = 0.0;
				float steps = g_max(1.0, distance / resolution);
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
					for(size_t i = 0; i < (size_t) rows * cols; ++i) 
						if(fillGrid[i] == 2) fillGrid.set(i, 0);
					step += 1.0;
				} while(found && step <= steps);
			}

			/**
			 * Blends two rasters together using the alpha grid for blending.
			 */
			void blend(Grid<float> &imgGrid, Grid<float> &bgGrid, Grid<float> &alpha, 
				int cols, int rows, float imNodata, float bgNodata) {
				for(int r = 0; r < rows; ++r) {
					for(int c = 0; c < cols; ++c) {
						float bv = bgGrid.get(c, r);
						float iv = imgGrid.get(c, r);
						if(!(bv == bgNodata)) {
							float av = alpha.get(c, r);
							if(iv == imNodata)
								iv = 0.0;
							bgGrid.set(c, r, bv * (1.0 - av) + iv * av);
						}
					}
				}
			}

		} // util

		/**
		 * Mosaic the given files together using the first as the base. The base file will serve as a clipping
		 * mask for the others. The spatial reference system and resolution of all layers must match.
		 * The distance argument determines the distance over which the edges are feathered. The
		 * overviews argument, if true, forces the construction of new overviews. If this is not
		* used, the existing overviews may obscure the fact that the image has changed (when zoomed out.)
		 */
		 // TODO: Background must be larger than other files. Remedy this.
		void mosaic(std::vector<std::string> &files, std::string &outfile, float distance, int rowHeight = 500, int threads = 1) {
			
			g_debug("mosaic");

			if(distance <= 0.0)
				g_argerr("Invalid distance: " << distance);
			if(outfile.size() == 0)
				g_argerr("No output file given.");
			if(files.size() < 2)
				g_argerr("Less than 2 files. Nothing to do.");

			using namespace geotools::raster::util;

			// Open the BG file for reading only.
			g_debug(" -- opening base file.");
			Raster<float> base(files[0]);

			// Create the destination file to modify it.
			Raster<float> output(outfile, 1, base);
			{
				MemRaster<float> buf(base.cols(), 1000);
				g_debug(" -- writing base file to output.");
				// Copy by block to avoid memory problems on big rasters. 
				// TODO: Configurable block size.
				for(int r = 0; r < base.rows(); r += 1000) {
					//base.readBlock(0, r, buf);
					//output.writeBlock(0, r, buf);
				}
			}

			// Determine if the vertical and horizontal resolution are positive.
			// Get the smaller resolution to use for distance calculations.
			bool posYRes = base.resolutionY() > 0;
			bool posXRes = base.resolutionX() > 0;
			float outNodata = base.nodata();
			float res = g_min(g_abs(base.resolutionX()), g_abs(base.resolutionY()));

			// Snap the distance to resolution
			distance = std::floor(distance / res) * res;
			g_debug(" -- snapped distance is " << distance);

			// The number of rows required to accomodate the fade distance.
			int buffer = (int) distance / res + 1;
			while(rowHeight <= buffer) {
				g_warn("The row height (" << rowHeight << ") must be larger than the buffer distance (" << buffer << "). Doubling.");
				rowHeight *= 2;
			}

			// Iterate over the files, adding each one to the background.
			for(unsigned int i = 1; i < files.size(); ++i) {

				g_debug(" -- processing file: " << files[i]);
				Raster<float> input(files[i]);
				if(output.resolutionX() != input.resolutionX() || output.resolutionY() != input.resolutionY()) 
					g_argerr("Raster's resolution does not match the background.");

				// Get the origin of the input w/r/t the output and the number of overlapping cols, rows.
				int col = output.toCol(posXRes ? input.minx() : input.maxx());
				int row = output.toRow(posYRes ? input.miny() : input.maxy());
				int cols = output.toCol(posXRes ? input.maxx() : input.minx()) - col;
				int rows = output.toRow(posYRes ? input.maxy() : input.miny()) - row;
				float imNodata = input.nodata();

				g_debug(" -- col: " << col << ", row: " << row << ", cols: " << cols << ", rows: " << rows);

				const int batchSize = ((int) std::ceil((double) rows / threads / rowHeight)) * rowHeight;

				g_debug(" -- batch size: " << batchSize << " for " << threads << " threads");

				#pragma omp parallel for
				for(int o = 0; o < threads; ++o) {
					// TODO: This is to solve a problem with skipping the second iteration under OMP.
					int offset = o * batchSize;
					g_debug(" -- o: " << o << "; offset: " << offset << "; batchSize: " << batchSize);

					// Initialize the grids.
					MemRaster<float> imGrid(cols + buffer * 2, rowHeight + buffer * 2);
					MemRaster<float> outGrid(cols + buffer * 2, rowHeight + buffer * 2);
					MemRaster<float> alphaGrid(cols + buffer * 2, rowHeight + buffer * 2);

					int inCol   = (col - buffer) < 0 ? -(col - buffer) : 0;
					int baseCol = (col - buffer) < 0 ? 0               : (col - buffer);
	
					for(int bufRow = offset; bufRow < offset + batchSize; bufRow += rowHeight) {

						g_debug(" -- bufrow: " << bufRow << "; row height: " << rowHeight << "; thread: " << omp_get_thread_num());
						
						// Set to zero or nodata depending on the band.
						outGrid.fill(outNodata);
						imGrid.fill(imNodata);
						alphaGrid.fill(1.0);

						// Load the overlay.
						g_debug(" -- loading overlay; thread: " << omp_get_thread_num());
						#pragma omp critical(access_input)
						{
							int a = bufRow - buffer < 0 ? 0 : bufRow - buffer;
							int b = bufRow - buffer < 0 ? -(bufRow - buffer) : 0;
							input.readBlock(inCol, a, imGrid, buffer, b);
						}

						// Feather the overlay
						if(distance > 0.0) {
							g_debug(" -- feathering; thread: " << omp_get_thread_num());
							int b = g_min(rowHeight + buffer * 2, rows - bufRow);
							feather(imGrid, alphaGrid, cols + buffer * 2, b, distance, imNodata, res);
						}

						// Read background data.
						g_debug(" -- reading background; thread: " << omp_get_thread_num());
						#pragma omp critical(access_output)
						{
							int a = bufRow + row - buffer < 0 ? 0 : bufRow + row - buffer;
							int b = bufRow + row - buffer < 0 ? -(bufRow + row - buffer) : 0;
							output.readBlock(baseCol, a, outGrid, 0, b);
						}

						// Blend the overlay into the output.
						g_debug(" -- blending; thread: " << omp_get_thread_num());
						//blend(imGrid, outGrid, alphaGrid, cols + buffer * 2, rowHeight + buffer * 2, imNodata, outNodata);

						// Write back to the output.
						// We are extracting a slice out of the buffer, not writing the whole thing.
						g_debug(" -- writing output; thread: " << omp_get_thread_num());
						#pragma omp critical(access_output)
						{
							int b = g_min(rowHeight, output.rows() - bufRow + row);
							output.writeBlock(baseCol, bufRow + row, outGrid, 0, buffer, 0, b);
						}
					}
				}
			}
		}

	} // raster

} // geotools

void usage() {
	std::cerr << "Usage: mosaic [options] -o <output file> <file [file [file [...]]]>\n"
		<< "    -o <file>        The output file.\n"
		<< "    -d <distance>    The feather distance in map units (default 100.)\n"
		<< "    <file [...]>     A list of files. The first is used as the background\n"
		<< "                     and determines the size of the output. Subsequent\n"
		<< "                     images are layered on top and feathered.\n"
		<< "    -v               Verbose messages.\n"
		<< "    -t               The number of threads to use. Defaults to the number\n"
		<< "                     of cores.\n"
		<< "    -r               Buffer row height. This manages the heights of the buffer\n"
		<< "                     that is used for feathering/blending. Keep in mind the width\n"
		<< "                     of the images, the number of threads and the size of the type.\n";
	
}

int main(int argc, char **argv) {

 	try {

	 	float distance = 100.0;
	 	std::vector<std::string> files;
	 	std::string outfile;
	 	int threads = 0;
	 	int rowHeight = 500;

	 	for(int i = 1; i < argc; ++i) {
	 		std::string arg(argv[i]);
	 		if(arg == "-d") {
	 			distance = atof(argv[++i]);
	 		} else if(arg == "-o") {
	 			outfile = argv[++i];
			} else if(arg == "-v") {
				g_loglevel(G_LOG_DEBUG);
			} else if(arg == "-t") {
				threads = atoi(argv[++i]);
			} else if(arg == "-r") {
				rowHeight = atoi(argv[++i]);
	 		} else {
	 			files.push_back(argv[i]);
	 		}
	 	}

	 	if(threads > 0) {
	 		g_debug("running with " << threads << " threads");
	 		omp_set_dynamic(1);
	 		omp_set_num_threads(threads);
	 	}

 		geotools::raster::mosaic(files, outfile, distance, rowHeight, threads);

 	} catch(const std::exception &e) {
 		g_error(e.what());
 		usage();
 		return 1;
 	}

 	return 0;
 }
