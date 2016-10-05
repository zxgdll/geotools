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
#include "mosaic.hpp"
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
			void blend(Grid<float> &imGrid, Grid<float> &bgGrid, Grid<float> &alpha, float imNodata, float bgNodata, int buffer) {
				for(int r = buffer; r < imGrid.rows() - buffer; ++r) {
					for(int c = buffer; c < imGrid.cols() - buffer; ++c) {
						float bv = bgGrid.get(c, r);
						float iv = imGrid.get(c, r);
						if(!(bv == bgNodata || iv == imNodata)) {
							float av = alpha.get(c, r);
							bgGrid.set(c, r, bv * (1.0 - av) + iv * av);
						}
					}
				}
			}

			int __tile_id = 0;

			class Tile {
			public:
				int id;
				int tileSize;
				int buffer;
				int iCol;
				int iRow;
				int oCol;
				int oRow;

				Tile(int tileSize, int buffer, int iCol, int iRow, int oCol, int oRow) :
					id(++__tile_id),
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

		} // util


		void Mosaic::setFileCallback(void (*callback)(float)) {
			m_fileCallback = callback;
		}

		void Mosaic::setOverallCallback(void (*callback)(float)) {
			m_overallCallback = callback;
		}

		/**
		 * Mosaic the given files together using the first as the base. The base file will serve as a clipping
		 * mask for the others. The spatial reference system and resolution of all layers must match.
		 * The distance argument determines the distance over which the edges are feathered. The
		 * overviews argument, if true, forces the construction of new overviews. If this is not
		* used, the existing overviews may obscure the fact that the image has changed (when zoomed out.)
		 */
		 // TODO: Background must be larger than other files. Remedy this.
		void Mosaic::mosaic(const std::vector<std::string> &files, const std::string &outfile, float distance, int tileSize, int threads) {
			
			g_debug("mosaic");

			if(distance <= 0.0)
				g_argerr("Invalid distance: " << distance);
			if(outfile.size() == 0)
				g_argerr("No output file given.");
			if(files.size() < 2)
				g_argerr("Less than 2 files. Nothing to do.");
			if(tileSize <= 0)
				g_argerr("Tile size must be greater than zero.");

		 	if(threads > 0) {
		 		g_debug("running with " << threads << " threads");
		 		omp_set_dynamic(1);
		 		omp_set_num_threads(threads);
		 	} else {
		 		g_argerr("Run with >=1 thread.");
		 	}

			using namespace geotools::raster::util;

			// Open the BG file for reading only.
			g_debug(" -- opening base file.");
			Raster<float> base(files[0]);
	
			// Check some properties of other files for compatibility.
			for(unsigned int i = 1; i < files.size(); ++i) {
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

			for(unsigned int i = 1; i < files.size(); ++i) {

				if(m_overallCallback)
					m_overallCallback((i - 0.5) / (files.size() - 1));

				Raster<float> input(files[i]);

				Bounds bounds = base.bounds().intersection(input.bounds());

				// Compute the column/row bounds of the intersection.
				int iStartCol = input.toCol(input.positiveX() ? bounds.minx() : bounds.maxx());
				int iStartRow = input.toRow(input.positiveY() ? bounds.miny() : bounds.maxy());
				int iEndCol = input.toCol(input.positiveX() ? bounds.maxx() : bounds.minx());
				int iEndRow = input.toRow(input.positiveY() ? bounds.maxy() : bounds.miny());
				int oStartCol = base.toCol(output.positiveX() ? bounds.minx() : bounds.maxx());
				int oStartRow = base.toRow(output.positiveY() ? bounds.miny() : bounds.maxy());

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

				std::vector<std::unique_ptr<Tile> > tiles;

				for(int ir = iStartRow, rr = oStartRow; ir <= iEndRow; ir += tileSize, rr += tileSize) {
					for(int ic = iStartCol, oc = oStartCol; ic <= iEndCol; ic += tileSize, oc += tileSize) {
						std::unique_ptr<Tile> t(new Tile(tileSize, buffer, ic, ir, oc, rr));
						tiles.push_back(std::move(t));
					}
				}

				int tileStatus = 0;

				#pragma omp parallel shared(tileStatus)
				{
					// Initialize the grids.
					MemRaster<float> inGrid(tileSize + buffer * 2, tileSize + buffer * 2);
					MemRaster<float> outGrid(tileSize + buffer * 2, tileSize + buffer * 2);
					MemRaster<float> alphaGrid(tileSize + buffer * 2, tileSize + buffer * 2);
					float outNodata = output.nodata();
					float inNodata = input.nodata();
					
					#pragma omp for nowait
					for(unsigned int t = 0; t < tiles.size(); ++t) {

						#pragma omp atomic
						tileStatus++;

						if(m_fileCallback)
							m_fileCallback((tileStatus - 0.5) / tiles.size());

						std::unique_ptr<Tile> tile = std::move(tiles[t]);

						g_debug(" -- tile " << t << ", " << tile->id << "; thread: " << omp_get_thread_num());

						// Set to zero or nodata depending on the band.
						outGrid.fill(outNodata);
						inGrid.fill(inNodata);
						alphaGrid.fill(1.0);

						// Load foreground
						if(!tile->readInput(inGrid, input))
							continue;

						// Feather returns true if it did any work.
						if(!feather(inGrid, alphaGrid, distance, inNodata, res)) 
							continue;

						// Load background
						if(!tile->readOutput(outGrid, output))
							continue;

						// Blend the fore/background images with alpha
						blend(inGrid, outGrid, alphaGrid, inNodata, outNodata, buffer);
						
						// Write to output.
						tile->writeOutput(outGrid, output);

						if(m_fileCallback)
							m_fileCallback((float) tileStatus / tiles.size());

					}

				} // parallel

				if(m_fileCallback)
					m_fileCallback(1.0);
				if(m_overallCallback)
					m_overallCallback((float) i / (files.size() - 1));
			}

			g_debug(" -- done.");

		}

	} // raster

} // geotools

