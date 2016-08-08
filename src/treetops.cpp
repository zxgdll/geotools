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

#include "trees.hpp"

using namespace geotools::raster;

namespace trees {

	namespace util {

		/**
		 * A simple class for maintaining information about a tree top.
		 */
		class Top {
		public:
			size_t m_id;
			double m_x;
			double m_y;
			double m_z;
			Top(size_t id, double x, double y, double z) :
				m_id(id),
				m_x(x), m_y(y), m_z(z) {
			}
			Top(const Top &t) :
				Top(t.m_id, t.m_x, t.m_y, t.m_z) {
			}
			Top() :
				Top(0, 0, 0, 0) {
			}
		};

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
	void treetops(std::string &inraster, std::string &topshp, int window) {

		if(inraster.empty())
			_argerr("Input raster cannot be empty.");
		if(topshp.empty())
			_argerr("The treetop output name cannot be empty.");
		if(window < 3)
			_argerr("A window size less than 3 makes no sense.");
		if(window % 2 == 0) {
			_warn("Window is " << window << ". Bumping up to " << (window + 1));
			window++;
		}

		using namespace trees::util;

		// Input raster.
		std::map<size_t, std::unique_ptr<Top> > tops;

		Raster<float> raster(inraster);
		size_t n = 0;
		size_t count = raster.size();
		int cols = raster.cols();
		int rows = raster.rows();
		int offset = window / 2;
		int row = 0;

		#pragma omp parallel
		{
			_trace("Thread: " << omp_get_thread_num() << "/" << omp_get_num_threads());
			MemRaster<float> blk(cols, window);
			float max;
			std::map<size_t, std::unique_ptr<Top> > tops0;
		
			while(row < rows) {
				for(int r = row; r < _min(row + window, rows - window); ++r) { 
					blk.nodata(raster.nodata());
					#pragma omp critical
					{
						raster.readBlock(0, r, blk);
					}
					for(int c = 0; c < cols - window; ++c) {
						size_t id = ((size_t) c << 32) | r;
						if(isMaxCenter(blk, c, 0, window, &max)) {
							std::unique_ptr<Top> t(new Top(id, 
								raster.toX(c + offset) + raster.resolutionX() / 2.0, // center of pixel
								raster.toY(r + offset) + raster.resolutionY() / 2.0, 
								max
							));
							tops0[id] = std::move(t);
						}
					}
				}

				#pragma omp atomic
				row += window;

				_status(row, rows);
			}

			_trace("Copying tops to output dict.");
			#pragma omp critical
			{
				for(auto it = tops0.begin(); it != tops0.end(); ++it) 
					tops[it->first] = std::move(it->second);
			}
		}

		_status(rows, rows, true);

		std::stringstream out;
		out << std::setprecision(12);
		out << "id,x,y,z" << std::endl;
		for(auto it = tops.begin(); it != tops.end(); ++it) {
			Top *t = it->second.get();
			out << t->m_id << "," << t->m_x << "," << t->m_y << "," << t->m_z << std::endl;
		}
		std::cout << out.str();
	}

}
