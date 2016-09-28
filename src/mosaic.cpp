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
			bool feather(Grid<float> &srcGrid, Grid<float> &dstGrid, float distance, float nodata, float resolution) {
				// Fill grid is used to keep track of where the edges are as they're "snowed in"
				// Starts out as a mask of non-nodata pixels from the source.
				int cols = srcGrid.cols();
				int rows = srcGrid.rows();
				//g_debug("feather: nodata: " << nodata << "; res: " << resolution << "; distance: " << distance << "; cols: " << cols << "; rows: " << rows);
				MemRaster<char> fillGrid(cols, rows);
				int valid = 0;
				for(size_t i = 0; i < (size_t) rows * cols; ++i) {
					if(srcGrid[i] == nodata) {
						fillGrid.set(i, 0);
					} else {
						fillGrid.set(i, 1);
						++valid;
					}
				}
			
				if(valid == 0) 
					return false;

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

				return true;
			}

			/**
			 * Blends two rasters together using the alpha grid for blending.
			 */
			void blend(Grid<float> &imGrid, Grid<float> &bgGrid, Grid<float> &alpha, float imNodata, float bgNodata) {
				for(int r = 0; r < imGrid.rows(); ++r) {
					for(int c = 0; c < imGrid.cols(); ++c) {
						float bv = bgGrid.get(c, r);
						float iv = imGrid.get(c, r);
						if(!(bv == bgNodata)) {
							float av = alpha.get(c, r);
							if(iv != imNodata)
								bgGrid.set(c, r, bv * (1.0 - av) + iv * av);
						}
					}
				}
			}

		} // util

		class Tile {
		public:
			int iCol;
			int iRow;
			int oCol;
			int oRow;
			int tileSize;
			int buffer;

			Tile(int tileSize, int buffer, int iCol, int iRow, int oCol, int oRow) :
				tileSize(tileSize), buffer(buffer),
				iCol(iCol), iRow(iRow),
				oCol(oCol), oRow(oRow) {
			}

			bool readInput(MemRaster<float> &buf, Raster<float> &input) const {
				int col = iCol - buffer;
				int row = iRow - buffer;
				int cols = tileSize + buffer * 2;
				int rows = tileSize + buffer * 2;
				int cOff = 0;
				int rOff = 0;
				if(col < 0) {
					cOff = -col;
					cols -= col;
					col = 0;
				}
				if(row < 0) {
					rOff = -row;
					rows -= row;
					row = 0;
				}
				if(cols <= 0 || rows <= 0)
					return false;
				input.readBlock(col, row, buf, cOff, rOff, cols, rows);
				return true;
			}
				
			bool readOutput(MemRaster<float> &buf, Raster<float> &output) const {
				int col = oCol - buffer;
				int row = oRow - buffer;
				int cols = tileSize + buffer * 2;
				int rows = tileSize + buffer * 2;
				int cOff = 0;
				int rOff = 0;
				if(col < 0) {
					cOff = -col;
					cols -= col;
					col = 0;
				}
				if(row < 0) {
					rOff = -row;
					rows -= row;
					row = 0;
				}
				if(cols <= 0 || rows <= 0)
					return false;
				output.readBlock(col, row, buf, cOff, rOff, cols, rows);
				return true;
			}

			void writeOutput(MemRaster<float> &buf, Raster<float> &output) const {
				if(!(oCol >= output.cols() || oRow >= output.rows()))
					output.writeBlock(oCol, oRow, buf, buffer, buffer, tileSize, tileSize);
			}

			void print() const {
				std::cerr << "[Tile: in: " << iCol << "," << iRow << "; out: " << oCol << "," << oRow << "]" << std::endl;
			}
		};

		/**
		 * Mosaic the given files together using the first as the base. The base file will serve as a clipping
		 * mask for the others. The spatial reference system and resolution of all layers must match.
		 * The distance argument determines the distance over which the edges are feathered. The
		 * overviews argument, if true, forces the construction of new overviews. If this is not
		* used, the existing overviews may obscure the fact that the image has changed (when zoomed out.)
		 */
		 // TODO: Background must be larger than other files. Remedy this.
		void mosaic(std::vector<std::string> &files, std::string &outfile, float distance, int tileSize, int threads = 1) {
			
			g_debug("mosaic");

			if(distance <= 0.0)
				g_argerr("Invalid distance: " << distance);
			if(outfile.size() == 0)
				g_argerr("No output file given.");
			if(files.size() < 2)
				g_argerr("Less than 2 files. Nothing to do.");
			if(tileSize <= 0)
				g_argerr("Tile size must be greater than zero.");

			using namespace geotools::raster::util;

			// Open the BG file for reading only.
			g_debug(" -- opening base file.");
			Raster<float> base(files[0]);
	
			// Check some properties of other files for compatibility.
			for(int i = 1; i < files.size(); ++i) {
				Raster<float> check(files[i]);
				if(check.resolutionX() != base.resolutionX() || check.resolutionY() != base.resolutionY())
					g_argerr("Resolution of " << files[i] << " doesn't match base.")
			}
				
			// Create the destination file to modify it.
			Raster<float> output(outfile, 1, base);
			{
				MemRaster<float> buf(base.cols(), 1000);
				g_debug(" -- writing base file to output.");
				// Copy by block to avoid memory problems on big rasters. 
				// TODO: Configurable block size.
				for(int r = 0; r < base.rows(); r += 1000) {
					base.readBlock(0, r, buf);
					output.writeBlock(0, r, buf);
				}
			}

			for(int i = 1; i < files.size(); ++i) {

				Raster<float> input(files[i]);

				Bounds bounds = base.bounds().intersection(input.bounds());

				// Compute the column/row bounds of the intersection.
				int iStartCol = input.toCol(input.positiveX() ? bounds.minx() : bounds.maxx());
				int iStartRow = input.toRow(input.positiveY() ? bounds.miny() : bounds.maxy());
				int iEndCol = input.toCol(input.positiveX() ? bounds.maxx() : bounds.minx());
				int iEndRow = input.toRow(input.positiveY() ? bounds.maxy() : bounds.miny());
				int oStartCol = base.toCol(output.positiveX() ? bounds.minx() : bounds.maxx());
				int oStartRow = base.toRow(output.positiveY() ? bounds.miny() : bounds.maxy());
				int oEndCol = base.toCol(output.positiveX() ? bounds.maxx() : bounds.minx());
				int oEndRow = base.toRow(output.positiveY() ? bounds.maxy() : bounds.miny());

				// Get the minimum absolute resolution.			
				double res = g_min(g_abs(base.resolutionX()), g_abs(base.resolutionY()));
				g_debug(" -- resolution is: " << res);

				// Snap the distance to resolution
				distance = std::floor(distance / res) * res;
				g_debug(" -- snapped distance is " << distance);

				// The number of rows required to accomodate the fade distance.
				int buffer = (int) distance / res + 1;
				while(tileSize <= buffer) {
					g_warn("The tileSize (" << tileSize << ") must be larger than the buffer distance (" << buffer << "). Doubling.");
					tileSize *= 2;
				}
				g_debug(" -- tile size is: " << tileSize);

				std::vector<Tile> tiles;

				for(int ir = iStartRow, rr = oStartRow; ir <= iEndRow; ir += tileSize, rr += tileSize) {
					for(int ic = iStartCol, oc = oStartCol; ic <= iEndCol; ic += tileSize, oc += tileSize)
						tiles.push_back(Tile(tileSize, buffer, ic, ir, oc, rr));
				}

				int t = 0;

				#pragma omp parallel
				{

					// Initialize the grids.
					MemRaster<float> inGrid(tileSize + buffer * 2, tileSize + buffer * 2);
					MemRaster<float> outGrid(tileSize + buffer * 2, tileSize + buffer * 2);
					MemRaster<float> alphaGrid(tileSize + buffer * 2, tileSize + buffer * 2);
					float outNodata = output.nodata();
					float inNodata = input.nodata();

					#pragma omp for
					for(int t = 0; t < tiles.size(); ++t) {
						const Tile &tile = tiles[t];
						g_debug(" -- tile " << t << "; thread: " << omp_get_thread_num());
		
						// Set to zero or nodata depending on the band.
						outGrid.fill(outNodata);
						inGrid.fill(inNodata);
						alphaGrid.fill(1.0);

						bool good;

						#pragma omp critical(input)
						{
							good = tile.readInput(inGrid, input);
						}
					
						if(!good) continue;

						#pragma omp critical(output)
						{
							good = tile.readOutput(outGrid, output);
						}

						if(!good) continue;

						// Feather returns true if it did any work.
						if(feather(inGrid, alphaGrid, distance, inNodata, res)) {
							blend(inGrid, outGrid, alphaGrid, inNodata, outNodata);
							#pragma omp critical(output)
							{
								tile.writeOutput(outGrid, output);
							}
						}
					}

				} // parallel
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
		<< "    -s               Tile size. This manages the sizes of the tiles. Default 1024.\n"
		<< "                     that is used for feathering/blending. Keep in mind the width\n"
		<< "                     of the images, the number of threads and the size of the type.\n";
	
}

int main(int argc, char **argv) {

 	try {

	 	float distance = 100.0;
	 	std::vector<std::string> files;
	 	std::string outfile;
	 	int threads = 0;
	 	int tileSize = 1024;

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
			} else if(arg == "-s") {
				tileSize = atoi(argv[++i]);
	 		} else {
	 			files.push_back(argv[i]);
	 		}
	 	}

	 	if(threads > 0) {
	 		g_debug("running with " << threads << " threads");
	 		omp_set_dynamic(1);
	 		omp_set_num_threads(threads);
	 	}

 		geotools::raster::mosaic(files, outfile, distance, tileSize, threads);

 	} catch(const std::exception &e) {
 		g_error(e.what());
 		usage();
 		return 1;
 	}

 	return 0;
 }
