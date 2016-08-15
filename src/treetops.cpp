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
#include "Util.hpp"
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

	}

	/**
	*/
	void treetops(std::string &inraster, std::string &outvect, std::map<size_t, std::unique_ptr<trees::util::Top> > &tops, int window) {

		if(inraster.empty())
			g_argerr("Input raster cannot be empty.");
		if(outvect.empty())
			g_warn("The treetop output filename is empty; not writing treetops");
		if(window < 3)
			g_argerr("A window size less than 3 makes no sense.");
		if(window % 2 == 0) {
			g_warn("Window is " << window << ". Bumping up to " << (window + 1));
			window++;
		}

		using namespace trees::util;

		Raster<float> raster(inraster);
		size_t n = 0;
		size_t count = raster.size();
		int cols = raster.cols();
		int rows = raster.rows();
		int offset = window / 2;
		int row = 0;

		SQLite db(outvect, SQLite::POINT, 26910, {{"id", 1}});
		
		// TODO: Set up output raster.

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
						if(isMaxCenter(blk, c, r, window, &max)) {
							std::unique_ptr<Top> t(new Top(id, 
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
		int c, r;
		Node(int c, int r) {
			this->c;
			this->r;
		}
	};

	void delineateCrown(Raster<float> &inrast, Raster<unsigned int> &outrast, trees::util::Top *top, double threshold) {
	
		std::queue<Node> q;
		std::vector<bool> visited((size_t) inrast.cols() * inrast.rows());

		q.push(Node(top->m_col, top->m_row));

		while(q.size()) {

			Node n = q.front();
			q.pop();

			visited[(size_t) n.r * inrast.cols() + n.c] = true;
			
			for(int r = g_max(0, n.r - 1); r < g_min(inrast.rows(), n.r + 1); ++r) {
			for(int c = g_max(0, n.c - 1); c < g_min(inrast.cols(), n.c + 1); ++c) {
		
				size_t idx = (size_t) r * inrast.cols() + c;
				if(c != top->m_col && r != top->m_row && !visited[idx]) {
				
					double v = inrast.get(c, r);
					double dif = v - inrast.get(top->m_col, top->m_row);
					if(dif < 0 && g_abs(dif / top->m_z) <= threshold) {
						visited[idx] = true;
						q.push(Node(c, r));
						outrast.set(c, r, top->m_id);
					}

				}
			}
			}
		}

	}

	void treecrowns(const std::string &infile, const std::string &outrfile, const std::string &outvfile, 
		std::map<size_t, std::unique_ptr<trees::util::Top> > &tops, double threshold) {

		Raster<float> inrast(infile);
		Raster<unsigned int> outrast(outrfile, inrast);
		outrast.nodata(0);
		outrast.fill(0);

		size_t n = 0;
		for(auto it = tops.begin(); it != tops.end(); ++it) {

			delineateCrown(inrast, outrast, it->second.get(), threshold);
			Util::status(n++, tops.size(), false);
		}

		Util::status(tops.size(), tops.size(), true);

	}


}
