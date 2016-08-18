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

namespace trees {

	namespace util {

		/**
		 * Returns true if the pixel at the center of the given
		 * raster is the maximum value in the raster.
		 */
		bool isMaxCenter(MemRaster<float> &raster, int col, int row, int window, float *max) {
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

	} // Util

	/**
	 * Locates tree top points on a canopy height model.
	 * inraster   - The input raster. If the smoothed parameter is given, it is assumed that the input
	 * 		raster is not smoothed, and smoothing is performed. Otherwise the raster is used 
	 * 		as-is.
	 * outvect    - The name of a vector file (spatialite) to save the tops to. The table name is "data"
	 * 		and the columns are "geom" and "id." Geom is a 3D point and ID is the top's unique ID.
	 * window     - The size of the kernel to use to locate tops. Greater to or equal to 3. Even values
	 * 		will be rounded up one.
	 * minHeight - The algorithm will not consider pixels below this height.
	 * smoothed   - Optional. A filename for the smoothed raster. If not given, inraster is used unchanged. 
	 * sigma      - Optional. If smoothed is given, the standard deviation for the gaussian kernel.
	 * kernel     - Optional. If smoothed is given, the kernel size for smoothing. 
	*/
	void treetops(const std::string &inraster, const std::string &outvect, std::map<size_t, std::unique_ptr<trees::util::Top> > &tops, 
		int window, double minHeight, const std::string &smoothed, double sigma, int kernel) {

		if(inraster.empty())
			g_argerr("Input raster cannot be empty.");
		if(outvect.empty())
			g_warn("The treetop output filename is empty; not writing treetops");
		if(window < 3)
			g_argerr("A window size less than 3 makes no sense.");
		if(!smoothed.empty()) {
			if(sigma <= 0)
				g_argerr("Sigma must be >0 for smoothing.");
			if(kernel < 3)
				g_argerr("Smoothing kernel must be 3 or larger.");
			if(kernel % 2 == 0) {
				++kernel;
				g_warn("Smoothing kernel must be an odd integer. Bumping up to " << kernel);
			}
		}
		if(window % 2 == 0) {
			window++;
			g_warn("Window is " << window << ". Bumping up to " << window);
		}

		using namespace trees::util;

		std::string rastfile;
		if(!smoothed.empty()) {
			g_trace("Smoothing raster: " << inraster << " -> " << smoothed);
			Raster<float> raster(inraster);
			Raster<float> sraster(smoothed, 1, raster);
			raster.smooth(sraster, sigma, kernel);
			rastfile.assign(smoothed);
		} else {
			g_warn("Using un-smoothed raster.");
			rastfile.assign(inraster);
		}

		Raster<float> raster(rastfile);

		size_t n = 0;
		size_t count = raster.size();
		int cols = raster.cols();
		int rows = raster.rows();
		int offset = window / 2;
		int row = 0;	
		size_t tid = 0;

		SQLite db(outvect, SQLite::POINT, 26910, {{"id", 1}});
		
		// This is the size of the cache used by each thread.
		int cachedRows = 500;
		int blockHeight = window + cachedRows;
		int topCount = 0;

		#pragma omp parallel
		{
			g_trace("Thread: " << omp_get_thread_num() << "/" << omp_get_num_threads());
			MemRaster<float> blk(cols, blockHeight);
			std::map<size_t, std::unique_ptr<Top> > tops0;
			float max;
			int topCount0 = 0;
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

				Util::status(g_min(row, rows * 2), rows * 2, false);
	
				#pragma omp atomic
				topCount += topCount0;
			}

			g_trace("Writing " << topCount << " tops.");
			//std::list<std::unique_ptr<Geom> > geoms;
			int tc = topCount / 2;
			#pragma omp critical(b)
			{
				db.begin();
				for(auto it = tops0.begin(); it != tops0.end(); ++it) {
					Top *t = it->second.get();
					db.addPoint(t->m_x, t->m_y, t->m_z, {{"id", std::to_string(t->m_id)}});
					tops[it->first] = std::move(it->second);
					if(++tc % 1000000 == 0) {
						Util::status(tc, topCount, false);
						db.commit();
						db.begin();
					}
				}
				db.commit();	
			}
		}

		Util::status(1, 1, true);
		g_trace("Done.");
	}

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
			Node(top->m_id, top->m_col, top->m_row, top->m_z, top->m_col, top->m_row, top->m_z) {
		}
	};

	double dist(int tc, int tr, int c, int r, double resolution) {
		return std::sqrt(g_sq((double) (tc - c)) + g_sq((double) (tr - r))) * resolution;
	}

	void delineateCrowns(Raster<float> &inrast, Raster<unsigned int> &outrast, const std::map<size_t, 
		std::unique_ptr<trees::util::Top> > &tops, double threshold, double radius) {

		std::queue<std::unique_ptr<Node> > q;
		for(auto it = tops.begin(); it != tops.end(); ++it)
			q.push(std::unique_ptr<Node>(new Node(it->second.get())));

		std::vector<bool> visited((size_t) inrast.cols() * inrast.rows());
		double nodata = inrast.nodata();
		double resolution = inrast.resolutionX();

		while(q.size()) {

			std::unique_ptr<Node> n = std::move(q.front());
			q.pop();

			size_t idx = (size_t) n->r * inrast.cols() + n->c;
			outrast.set(idx, n->id);

			for(int r = g_max(0, n->r - 1); r < g_min(inrast.rows(), n->r + 2); ++r) {
				for(int c = g_max(0, n->c - 1); c < g_min(inrast.cols(), n->c + 2); ++c) {
			//std::list<std::pair<int, int> > sites;
			//sites.push_back(std::pair<int,int>(n->r - 1, n->c));
			//sites.push_back(std::pair<int,int>(n->r + 1, n->c));
			//sites.push_back(std::pair<int,int>(n->r, n->c - 1));
			//sites.push_back(std::pair<int,int>(n->r, n->c + 1));

			// TODO: Option for d4/d8

			//for(auto it = sites.begin(); it != sites.end(); ++it) {
			//	int r = it->first;
			//	int c = it->second;
				if(r < 0 || c < 0 || r >= inrast.rows() || c >= inrast.cols()) continue;

					idx = (size_t) r * inrast.cols() + c;
					if(visited[idx])
						continue;
					double v = inrast.get(c, r);
					if(v != nodata && v < n->z && /*(v / n->tz) >= threshold &&*/ dist(n->tc, n->tr, n->c, n->r, resolution) <= radius) {
						q.push(std::unique_ptr<Node>(new Node(n->id, c, r, v, n->tc, n->tr, n->tz)));
						outrast.set(idx, n->id);
						visited[idx] = true;
					}
				}

			}
		}


	}

	/**
	 * infile    - The input raster file. This should be the same one used to delineate the 
	 *             treetop points; if a smooth raster was used for tops, it should be used for crowns.
	 * outrfile  - The output raster file. This file is the primary output, so not specifying it
	 *	       saves no resources. The output is an integer raster containing the IDs of the 
	 *	       treetop points which seed the algorithm.
	 * outvfile  - A vector file representing the crowns; extracted from the same raster as outrfile.
	 *	       This is an sqlite file; shapefiles have a size restriction.
	 * tops      - This is a map where the keys are a unique ID and the values are instances of Top, 
	 * 	       which contains information about the tree tops.
	 * threshold - A pixel is excluded if its height is within a specified proportion of the tree top's
	 * 	       height. This is a value between 0 and 1.
	 * radius    - Tree crowns will be clipped to this radius.
	 */
	void treecrowns(const std::string &infile, const std::string &outrfile, const std::string &outvfile, 
		std::map<size_t, std::unique_ptr<trees::util::Top> > &tops, double threshold, double radius) {

		Raster<float> inrast(infile);
		Raster<unsigned int> outrast(outrfile, 1, inrast);
		outrast.nodata(0);
		outrast.fill(0);

		delineateCrowns(inrast, outrast, tops, threshold, radius);

		Util::status(tops.size(), tops.size(), true);

	}


}
