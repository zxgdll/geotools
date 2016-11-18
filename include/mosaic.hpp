#ifndef __MOSAIC_HPP__
#define __MOSAIC_HPP__

#include <vector>
#include <string>

#include "geotools.h"
#include "util.hpp"

namespace geotools {
	
	namespace raster {

		class DLL_EXPORT Mosaic {
		private:
			geotools::util::Callbacks *m_callbacks;

		public:
			Mosaic();
			void setCallbacks(geotools::util::Callbacks *callbacks);
			void mosaic(const std::vector<std::string> &files, const std::string &outfile, float distance, int tileSize, int threads = 1);
		};

	} // raster

} // geotools


#endif