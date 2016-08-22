/*
 * treecrowns.cpp
 *
 *  Created on: May 3, 2016
 *      Author: rob
 */

#include <queue>
#include <iostream>
#include <omp.h>

#include "Raster.hpp"
#include "util.hpp"
#include "trees.hpp"

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
				int c, r, tc, tr;
				double z, tz;
				size_t id;
				Node(size_t id, int c, int r, double z, int tc, int tr, double tz) :
					id(id), 
					c(c), r(r), z(z), 
					tc(tc), tr(tr), tz(tz) {
				}
				Node(const trees::util::Top *top) :
					Node(top->id, top->col, top->row, top->z, top->col, top->row, top->z) {
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
				int mc, mr;
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

Top::Top(size_t id, double x, double y, double z, int col, int row) :
	id(id),
        x(x), y(y), z(z),
        col(col), row(row) {
}

Top::Top(const Top &t) :
	Top(t.id, t.x, t.y, t.z, t.col, t.row) {
}

Top::Top() :
	Top(0, 0, 0, 0, 0, 0) {
}


void TreeUtil::smooth(const std::string &inraster, const std::string &outraster, double sigma, double window) {
	if(inraster.empty())
		g_argerr("Input raster must be given.");
	if(outraster.empty())
		g_argerr("Output raster must be given.");
	Raster<float> in(inraster);
	Raster<float> out(outraster, 1, in);
	in.smooth(out, sigma, window);
}

void TreeUtil::treetops(const std::string &inraster, const std::string &outvect, int window, double minHeight,
	int srid, std::vector<std::unique_ptr<trees::util::Top> > *tops) {
	
	Util::status(0, 1, "Locating tops...");
	
	if(inraster.empty())
		g_argerr("Input raster cannot be empty.");
	if(outvect.empty())
		g_warn("The treetop output filename is empty; not writing treetops");
	if(window < 3)
		g_argerr("A window size less than 3 makes no sense.");
	if(window % 2 == 0) {
		window++;
		g_warn("Window is " << window << ". Bumping up to " << window);
	}

	Raster<float> raster(inraster);
	size_t count = raster.size();
	int cols = raster.cols();
	int rows = raster.rows();

	int offset = window / 2;
	int row = 0;	
	size_t tid = 0;

	//  TODO: Get SRID from raster.
	SQLite db(outvect, SQLite::POINT, srid, {{"id", 1}});
	db.makeFast();
	db.dropGeomIndex();
	db.setCacheSize(1024 * 1024);
	db.clear();
	
	// This is the size of the cache used by each thread. TODO: Make configurable.
	int cachedRows = 500;
	int blockHeight = window + cachedRows;
	size_t tc, topCount = 0;

	#pragma omp parallel
	{
		g_trace("Thread: " << omp_get_thread_num() << "/" << omp_get_num_threads());
		MemRaster<float> blk(cols, blockHeight);
		std::map<size_t, std::unique_ptr<Top> > tops0;
		double max;
		size_t topCount0 = 0;
		int curRow;

		while(true) {
			#pragma omp critical(a)
			{
				curRow = row;
				row += cachedRows;
			}
			if(curRow >= rows)
				break;
			#pragma omp critical(c)
			{
				blk.nodata(raster.nodata());
				raster.readBlock(0, curRow, blk);
			}
			for(int r = 0; r < cachedRows - window; ++r) {
				for(int c = 0; c < cols - window; ++c) {
					size_t id = ((size_t) c << 32) | (r + curRow);
					if(blk.get(c, r) >= minHeight && isMaxCenter(blk, c, r, window, &max)) {
						std::unique_ptr<Top> t(new Top(++tid, 
							raster.toX(c + offset) + raster.resolutionX() / 2.0, // center of pixel
							raster.toY(r + curRow + offset) + raster.resolutionY() / 2.0, 
							max,
							c + offset,
							r + curRow + offset
						));
						tops0[id] = std::move(t);
						++topCount0;
					}
				}
			}

			Util::status(row / 2, rows, "Locating tops...");
	
			#pragma omp atomic
			topCount += topCount0;
		}

		g_trace("Writing " << topCount << " tops.");

		Util::status(1, 2, "Locating tops...");

		size_t tc0, batch = db.maxAddPointCount();

		std::vector<std::unique_ptr<Point> > points;
		auto it = tops0.begin();
		for(size_t b = 0; b < topCount0; b += batch) {
			for(size_t b0 = 0; b0 < batch && it != tops0.end(); ++it, ++b0) {
				Top *t = it->second.get();
				Point *pt = new Point(t->x, t->y, t->z, {{"id", std::to_string(t->id)}});
				points.push_back(std::unique_ptr<Point>(pt));
				// Only add the result if the output vector is defined.
				if(tops)
					tops->push_back(std::move(it->second));
				++tc0;
			}

			#pragma omp atomic
			tc += tc0;
			Util::status(tc, topCount, "Locating tops...");

			#pragma omp critical(b)
			{
				g_trace("inserting " << points.size() << " points");
				db.begin();	
				// Add to the file.
				db.addPoints(points);
				db.commit();
			}
			points.clear();
		}
	}

	db.makeSlow();
	db.setCacheSize(1024);
	db.createGeomIndex();

	Util::status(1, 1, "Locating tops... Done.", true);
	g_trace("Done.");

}

void TreeUtil::treecrowns(const std::string &inraster, const std::string &topsvect, 
	const std::string &crownsrast, const std::string &crownsvect, 
	double threshold, double radius, double minHeight, bool d8) {

	Util::status(0, 1, "Delineating crowns...");

	Raster<float> inrast(inraster);
	Raster<unsigned int> outrast(crownsrast, 1, inrast);
	outrast.nodata(0);
	outrast.fill(0);
	double nodata = inrast.nodata();
	double resolution = inrast.resolutionX();

	Util::status(0, 1, "Delineating crowns...");

	// Initialize the database.
	SQLite db(topsvect);
	size_t geomCount;
	db.getGeomCount(&geomCount);
	g_trace("Processing " << geomCount << " tree tops.");

	// The number of extra rows above and below the buffer.
	int bufRows = (int) std::ceil(g_abs(radius / inrast.resolutionY()));
	// The height of the row, not including disposable buffer.
	int rowStep = (int) g_abs(std::ceil(100000.0 / geomCount * inrast.rows()) / inrast.resolutionY());
	//int rowStep = g_abs((geomCount / 100000) / inrast.resolutionY());
	// The totoal height of the buffer
	int rowHeight = rowStep + bufRows * 2;

	// Build the list of offsets for D8 or D4 search.
	std::vector<std::pair<int, int> > offsets; // pairs of col, row
	if(d8) {
		offsets.assign({{-1, -1}, {0, -1}, {1, -1}, {-1, 0}, {0, 0}, {1, 0}, {-1, 1}, {0, 1}, {1, 1}});
	} else {
		offsets.assign({{0, -1}, {-1, 0}, {1, 0}, {1, 1}});
	}

	#pragma omp parallel for
	for(int row = 0; row < inrast.rows(); row += rowStep) {

		g_trace("Processing 100,000 of " << geomCount << " points.");
		g_trace(" - row " << row << " of " << inrast.rows());

		MemRaster<unsigned int> blk(inrast.cols(), rowHeight);
		blk.fill(0);

		// Load the tree tops for the strip.
		Bounds bounds(inrast.toX(0), inrast.toY(row - bufRows), inrast.toX(inrast.cols()), inrast.toY(row + rowStep + bufRows * 2));
		std::vector<std::unique_ptr<Point> > tops;
		#pragma omp critical(a)
		{
			db.getPoints(tops, bounds);
		}

		// Convert the Tops to Nodes.
		std::queue<std::unique_ptr<Node> > q;
		for(auto it = tops.begin(); it != tops.end(); ++it) {
			Point *t = it->get();
			int col = inrast.toCol(t->x);
			int row = inrast.toRow(t->y);
			int id = atoi(t->fields["id"].c_str());
			q.push(std::unique_ptr<Node>(new Node(id, col, row, t->z, col, row, t->z)));
		}

		// To keep track of visited cells.
		std::vector<bool> visited((size_t) inrast.cols() * inrast.rows());

		// Run through the queue.
		while(q.size()) {

			std::unique_ptr<Node> n = std::move(q.front());
			q.pop();

			if(n->c >= 0 && n->c < blk.cols() && (n->r - row) >= 0 && (n->r - row) < blk.rows())
				blk.set(n->c, n->r - row, n->id);

			for(std::pair<int, int> offset : offsets) {
				int c = n->c + offset.first;
				int r = n->r + offset.second;
				
				if(r < 0 || c < 0 || r >= inrast.rows() || c >= inrast.cols()) continue;

				size_t idx = (size_t) r * inrast.cols() + c;
				if(visited[idx])
					continue;

				double v = inrast.get(c, r);
				if(v != nodata 								// is not nodata
					&& v < n->z 							// is less than the neighbouring pixel
					&& v >= minHeight 						// is greater than the min height
					&& (v / n->tz) >= threshold 					// is greater than the threshold height
					&& dist(n->tc, n->tr, n->c, n->r, resolution) <= radius		// is within the radius
				) {
					q.push(std::unique_ptr<Node>(new Node(n->id, c, r, v, n->tc, n->tr, n->tz)));
					visited[idx] = true;
				}
			}
		}

		if(row == 0) {
			MemRaster<unsigned int> tmp(blk.cols(), rowStep + bufRows);
			blk.readBlock(0, 0, tmp);
			#pragma omp critical(b)
			{
				outrast.writeBlock(0, row, tmp);
			}
		} else {
			MemRaster<unsigned int> tmp(blk.cols(), rowStep);
			blk.readBlock(0, bufRows, tmp);
			#pragma omp critical(c)
			{
				outrast.writeBlock(0, row + bufRows, tmp);
			}
		}

		Util::status(g_min(row, inrast.rows()), inrast.rows(), "Delineating crowns...");

	}

	Util::status(1, 1, "Delineating crowns... Done", true);


}
