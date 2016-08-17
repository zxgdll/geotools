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

		/**
		 * Compute the table of Gaussian weights given the size of the table
		 * and the std. deviation.
		 */
		void gaussianWeights(double *weights, int size, double sigma) {
			if(size % 2 == 0) ++size;
			for(int r = 0; r < size; ++r) {
				for(int c = 0; c < size; ++c) {
					int x = size / 2 - c;
					int y = size / 2 - r;
					weights[r * size + c] = (1 / (2 * G_PI * sigma * sigma)) * pow(G_E, -((x * x + y * y) / (2.0 * sigma * sigma)));
					//g_trace("weights " << c << ", " << r << "; " << x << ", " << y << "; " << weights[r * size + c]);
				}
			}
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

		/**
		 * Smooth the raster and write the smoothed version to the output raster.
		 */
		// TODO: No accounting for nodata.
		// TODO: Move to Raster.
		void smooth(Grid<float> &raster, Grid<float> &smoothed, double sigma, int size) {
			double weights[size * size];
			gaussianWeights(weights, size, sigma);
			float nodata = raster.nodata();
			for(int r = size / 2; r < raster.rows() - size / 2; ++r) {
				for(int c = size / 2; c < raster.cols() - size / 2; ++c) {
					double t = 0.0;
					double v;
					for(int gr = 0; gr < size; ++gr) {
						for(int gc = 0; gc < size; ++gc) {
							int rc = c - size / 2 + gc;
							int rr = r - size / 2 + gr;
							v = raster.get(rc, rr);
							if(v != nodata) {
								t += weights[gr * size + gc] * v;
							} else {
								g_runerr("Nodata found at " << rc << ", " << rr);
							}
						}
					}
					smoothed.set(c, r, t);
				}
			}
		}


	} // util

	/**
	*/
	void treetops(const std::string &inraster, const std::string &outvect, std::map<size_t, std::unique_ptr<trees::util::Top> > &tops, int window, const std::string &smoothed) {

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

		std::string rastfile;
		if(!smoothed.empty()) {
			g_trace("Smoothing raster: " << inraster << " -> " << smoothed);
			Raster<float> raster(inraster);
			Raster<float> sraster(smoothed, 1, raster);
			smooth(raster, sraster, 0.8, 3);
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
						if(isMaxCenter(blk, c, r, window, &max)) {
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
		int c, r;
		double z, tz;
		size_t id;
		Node(size_t _id, int _c, int _r, double _z, double _tz) :
			id(_id), 
			c(_c), 
			r(_r), 
			z(_z), 
			tz(_tz) {
		}
		Node(const trees::util::Top *top) :
			Node(top->m_id, top->m_col, top->m_row, top->m_z, top->m_z) {
		}
	};


	void delineateCrowns(Raster<float> &inrast, Raster<unsigned int> &outrast, const std::map<size_t, std::unique_ptr<trees::util::Top> > &tops, double threshold) {

		std::queue<std::unique_ptr<Node> > q;
		for(auto it = tops.begin(); it != tops.end(); ++it) {
			q.push(std::unique_ptr<Node>(new Node(it->second.get())));
			break;
		}

		std::vector<bool> visited((size_t) inrast.cols() * inrast.rows());
		size_t idx;
		double nodata = inrast.nodata();

		while(q.size()) {

			std::unique_ptr<Node> n = std::move(q.front());
			q.pop();

			g_trace("n " << n->c << ", " << n->r);
			idx = (size_t) n->r * inrast.cols() + n->c;

			outrast.set(idx, n->id);

			for(int r = g_max(0, n->r - 1); r < g_min(inrast.rows(), n->r + 2); ++r) {
				for(int c = g_max(0, n->c - 1); c < g_min(inrast.cols(), n->c + 2); ++c) {
					g_trace("node " << c << "," << r << "," << n->id);
					idx = (size_t) r * inrast.cols() + c;
					if(visited[idx])
						continue;
					double v = inrast.get(c, r);
					if(v != nodata && v < n->z) { // && (v / n->tz) <= threshold) {
						q.push(std::unique_ptr<Node>(new Node(n->id, c, r, v, n->tz)));
						outrast.set(idx, n->id);
						visited[idx] = true;
					}
				}
			}
			g_trace("q " << q.size());
		}


	}

	void treecrowns(const std::string &infile, const std::string &outrfile, const std::string &outvfile, 
		std::map<size_t, std::unique_ptr<trees::util::Top> > &tops, double threshold) {

		Raster<float> inrast(infile);
		Raster<unsigned int> outrast(outrfile, 1, inrast);
		outrast.nodata(0);
		outrast.fill(0);

		delineateCrowns(inrast, outrast, tops, threshold);

		Util::status(tops.size(), tops.size(), true);

	}


}
