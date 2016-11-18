/*
 * treecrowns.cpp
 *
 *  Created on: May 3, 2016
 *      Author: rob
 */

#include <queue>
#include <iostream>
#include <omp.h>

#include "raster.hpp"
#include "util.hpp"
#include "treetops.hpp"
#include "sqlite.hpp"

using namespace geotools::raster;
using namespace geotools::util;
using namespace geotools::db;

namespace geotools {

	namespace trees {

		namespace util {

			/**
			 * Represents a grid cell, and maintains some properties of the seed that originated it.
			 */
			class Node {
			public:
				uint64_t id;
				int c, r, tc, tr;
				double z, tz;
				Node(uint64_t id, int c, int r, double z, int tc, int tr, double tz) :
					id(id), 
					c(c), r(r), tc(tc), tr(tr),
					z(z), tz(tz) {
				}
			};

			/**
			 * Returns the distance between the given column/row pair in map units.
			 */
			double dist(int tc, int tr, int c, int r, double resolution) {
				return std::sqrt(g_sq((double) (tc - c)) + g_sq((double) (tr - r))) * resolution;
			}

			/**
			 * Returns true if the pixel at the center of the given
			 * raster is the maximum value in the raster.
			 */
			bool isMaxCenter(MemRaster<float> &raster, int col, int row, int window, double *max) {
				int cc = col + window / 2;
				int cr = row + window / 2;
				float nd = raster.nodata();
				if(raster.get(cc, cr)  == nd)
					return false;
				*max = 0;
				int mc = 0, mr = 0;
				for(int r = row; r < row + window; ++r) {
					for(int c = col; c < col + window; ++c) {
						float v = raster.get(c, r);
						if(v != nd && v > *max) {
							*max = v;
							mc = c;
							mr = r;
						}
					}
				}
				return mc == cc && mr == cr;
			}

			/**
			 * Guess a value for a cell, based on its neighbours.
			 */
			double interpNodata(Grid<float> &rast, int col, int row) {
				int size = 1;
				double nodata = rast.nodata();
				double v, t;
				int n;	
				while(size < 1000) {
					n = 0;
					t = 0;
					for(int c = g_max(0, col - size); c < g_min(rast.cols(), col + size + 1); ++c) {
						v = rast.get(c, g_max(0, row - size));
						if(v != nodata) t += v, ++n;
						v = rast.get(c, g_min(rast.rows() - 1, row + size));
						if(v != nodata) t += v, ++n;
					}
					for(int r = g_max(1, row - size); r < g_min(rast.rows(), row + size + 1); ++r) {
						v = rast.get(g_max(0, col - size), r);
						if(v != nodata) t += v, ++n;
						v = rast.get(g_min(rast.cols() - 1, col + size), r);
						if(v != nodata) t += v, ++n;
					}
					if(n > 0)
						return t / n;
					++size;
				}
				g_runerr("Couldn't find a pixel to use as fill.");
			}

		} // util

	} // trees

} // geotools

using namespace geotools::trees;
using namespace geotools::trees::util;
using namespace geotools::trees::config;

Top::Top(uint64_t id, double x, double y, double z, int col, int row) :
	id(id),
    x(x), y(y), z(z),
    col(col), row(row) {
}

Top::Top() :
	Top(0, 0, 0, 0, 0, 0) {
}


void Trees::smooth(const std::string &inraster, const std::string &outraster, double sigma, uint32_t window) {
	if(inraster.empty())
		g_argerr("Input raster must be given.");
	if(outraster.empty())
		g_argerr("Output raster must be given.");
	Raster<float> in(inraster);
	Raster<float> out(outraster, 1, in);
	in.smooth(out, sigma, window);
}

void Trees::treetops(const TreeTopConfig &config) {
	
	config.check();
	
	Util::status(0, 1, "Locating tops...");

	// Smooth if necessary; set appropriate input filename.
	std::string inputFilename;
	if(config.performSmoothing) {
		smooth(config.inputFilename, config.smoothedFilename, config.smoothingSigma, config.smoothingWindow);
		inputFilename = config.smoothedFilename;
	} else {
		inputFilename = config.inputFilename;
	}

	Util::status(0, 1, "Locating tops... [Preparing input raster.]");

	// Initialize input raster.
	Raster<float> raster(inputFilename);

	size_t row = 0;	
	size_t tid = 0;

	Util::status(0, 1, "Locating tops... [Preparing database.]");	

	// Prepare database.
	SQLite db(config.outputFilename, SQLite::POINT, config.srid, {{"id", 1}}); // TODO: Get SRID from raster.
	db.makeFast();
	g_debug("x");
	db.dropGeomIndex();
	db.setCacheSize(config.tableCacheSize);
	db.clear(); // TODO: Faster to delete it and start over.
	// This is the size of the cache used by each thread.
	int cachedRows = g_max(100, (config.rowCacheSize / raster.rows() / sizeof(float))); // TODO: Make row cache configurable.
	int blockHeight = config.searchWindow + cachedRows;
	uint64_t tc = 0;
	uint64_t topCount = 0;

	#pragma omp parallel
	{
		g_debug(" - thread: " << omp_get_thread_num() << "/" << omp_get_num_threads());
		MemRaster<float> blk(raster.cols(), blockHeight);
		std::map<uint64_t, std::unique_ptr<Top> > tops0;
		double max;
		uint64_t topCount0 = 0;
		int curRow;

		while(true) {
			
			Util::status(row / 2, raster.rows(), "Locating tops...");

			#pragma omp critical(a)
			{
				curRow = row;
				row += cachedRows;
			}
			
			if(curRow >= raster.rows())
				break;
			
			blk.nodata(raster.nodata());
			#pragma omp critical(c)
			{
				raster.readBlock(0, curRow, blk);
			}
			
			for(int r = 0; r < cachedRows - config.searchWindow; ++r) {
				for(int c = 0; c < raster.cols() - config.searchWindow; ++c) {
					uint64_t id = ((uint64_t) c << 32) | (r + curRow);
					if(blk.get(c, r) >= config.minHeight && isMaxCenter(blk, c, r, config.searchWindow, &max)) {
						std::unique_ptr<Top> t(new Top(
							++tid, 
							raster.toCentroidX(c + config.searchWindow / 2), // center of pixel
							raster.toCentroidY(r + curRow + config.searchWindow / 2), 
							max,
							c + config.searchWindow / 2,
							r + curRow + config.searchWindow / 2
						));
						tops0[id] = std::move(t);
						++topCount0;
					}
				}
			}
		}

		#pragma omp atomic
		topCount += topCount0; // TODO: Doesn't work for status because there's no barrier here.
		
		g_debug(" - writing " << topCount << " tops.");
		Util::status(1, 2, "Locating tops... [Saving to database.]");

		uint64_t b = 0;
		uint64_t batch = db.maxAddPointCount();

		std::vector<std::unique_ptr<Point> > points;
		for(auto it = tops0.begin(); it != tops0.end(); ++it) {
			
			const std::unique_ptr<Top> &t = it->second;
			
			std::unique_ptr<Point> pt(new Point(t->x, t->y, t->z, {{"id", std::to_string(t->id)}}));
			points.push_back(std::move(pt));

			if(++b % batch == 0 || b == topCount0) {	
				#pragma omp critical(b)
				{
					tc += points.size();
					g_debug("inserting " << points.size() << " points");
					db.begin();	
					db.addPoints(points);
					db.commit();
				}
				points.clear();
				Util::status(tc, topCount, "Locating tops... [Saving to database.]");
			}
		}
	}

	if(config.buildIndex) {
		Util::status(99, 100, "Locating tops... [Building index.]");
		db.createGeomIndex();
	}
	db.makeSlow();

	Util::status(1, 1, "Locating tops... Done.", true);
	g_debug("Done.");

}

void Trees::treecrowns(const std::string &inraster, const std::string &topsvect, 
	const std::string &crownsrast, const std::string &crownsvect, 
	double threshold, double radius, double minHeight, bool d8) {

	Util::status(0, 1, "Delineating crowns...");

	int batchSize = 1000000;

	Util::status(0, 1, "Delineating crowns... [Preparing input/output rasters.]");

	// Initialize the rasters.
	Raster<float> inrast(inraster);
	Raster<uint32_t> outrast(crownsrast, 1, inrast);
	outrast.nodata(0);
	outrast.fill(0);
	double nodata = inrast.nodata();
	double resolution = inrast.resolutionX();

	Util::status(0, 1, "Delineating crowns... [Preparing database.]");

	// Initialize the database, get the treetop cound and estimate the buffer size.
	SQLite db(topsvect);
	uint64_t geomCount;
	db.getGeomCount(&geomCount);
	g_debug("Processing " << geomCount << " tree tops.");

	// The number of extra rows above and below the buffer.
	int bufRows = (int) std::ceil(g_abs(radius / inrast.resolutionY()));
	// The height of the row, not including disposable buffer. Use bufRows as lower bound to
	// avoid read error later (keeps min row index to >=0)
	int rowStep = g_max(bufRows, (int) g_abs(std::ceil((double) batchSize / geomCount * inrast.rows()) / inrast.resolutionY()));
	// The totoal height of the buffer
	int rowHeight = rowStep + bufRows * 2;
	int rowCompleted = 0;

	// Build the list of offsets for D8 or D4 search.
	std::vector<std::pair<int, int> > offsets; // pairs of col, row
	if(d8) {
		offsets.assign({{-1, -1}, {0, -1}, {1, -1}, {-1, 0}, {0, 0}, {1, 0}, {-1, 1}, {0, 1}, {1, 1}});
	} else {
		offsets.assign({{0, -1}, {-1, 0}, {1, 0}, {1, 1}});
	}

	#pragma omp parallel for
	for(int row = 0; row < inrast.rows(); row += rowStep) {

		g_debug("Processing " << batchSize << " of " << geomCount << " points.");
		g_debug(" - row " << row << " of " << inrast.rows() << "; step: " << rowStep);

		Util::status(rowCompleted, inrast.rows(), "Delineating crowns...");
		
		// To keep track of visited cells.
		std::vector<bool> visited((uint64_t) inrast.cols() * rowHeight);
		MemRaster<float> buf(inrast.cols(), rowHeight);
		MemRaster<uint32_t> blk(inrast.cols(), rowHeight);
		buf.fill(inrast.nodata());
		blk.fill(0);

		// Load the tree tops for the strip.
		Bounds bounds(inrast.toX(0), inrast.toY(row - bufRows), inrast.toX(inrast.cols()), inrast.toY(row + rowStep + bufRows));
		std::vector<std::unique_ptr<Point> > tops;
		#pragma omp critical(crowns_getpoints)
		{
			db.getPoints(tops, bounds);
		}
		g_debug(" - reading from source: " << 0 << ", " << (row == 0 ? row : row - bufRows) << ", [buf], " << 0 << ", " << (row == 0 ? bufRows : 0));
		#pragma omp critical(crowns_readbuf) 
		{
			inrast.readBlock(0, row == 0 ? row : row - bufRows, buf, 0, row == 0 ? bufRows : 0);
		}

		// Convert the Tops to Nodes.
		std::queue<std::unique_ptr<Node> > q;
		for(const std::unique_ptr<Point> &t : tops) {
			int col = inrast.toCol(t->x);
			int row = inrast.toRow(t->y);
			int id = atoi(t->fields["id"].c_str());
			q.push(std::unique_ptr<Node>(new Node(id, col, row, t->z, col, row, t->z)));
		}

		// Run through the queue.
		while(q.size()) {

			std::unique_ptr<Node> n = std::move(q.front());
			q.pop();

			blk.set(n->c, n->r - row + bufRows, n->id);

			for(const std::pair<int, int> &offset : offsets) {
				int c = n->c + offset.first;
				int r = n->r + offset.second;
				
				//g_debug(" -- " << n->c << ", " << n->r << "; " << c << ", " << r);
				if(r < 0 || c < 0 || r >= inrast.rows() || c >= inrast.cols()) continue;
				if(r - row + bufRows < 0 || r - row  + bufRows >= buf.rows()) continue;

				uint64_t idx = (uint64_t) (r - row + bufRows) * inrast.cols() + c;
				if(visited[idx])
					continue;

				double v = buf.get(idx);
				if(v != nodata 													// is not nodata
					&& v < n->z 												// is less than the neighbouring pixel
					&& v >= minHeight 											// is greater than the min height
					&& (v / n->tz) >= threshold 								// is greater than the threshold height
					&& std::pow(n->tc - c, 2) + std::pow(n->tr - r, 2) <= std::pow(radius, 2)		// is within the radius
				) {
					q.push(std::unique_ptr<Node>(new Node(n->id, c, r, v, n->tc, n->tr, n->tz)));
					visited[idx] = true;
				}
			}
		}

		Util::status(rowCompleted, inrast.rows(), "Delineating crowns... [Writing to raster.]");
		g_debug(" - tmp block: cols: " << blk.cols() << ", rows: " <<  (row == 0 ? rowStep + bufRows : rowStep));
		if(row > 0 && (row + bufRows) >= inrast.rows())
			continue;
		#pragma omp critical(b)
		{
			outrast.writeBlock(0, row, blk, 0, bufRows);
			rowCompleted += rowStep;
		}
	}

	Util::status(1, 1, "Delineating crowns... Done", true);
}
