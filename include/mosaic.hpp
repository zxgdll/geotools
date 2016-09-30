#ifndef __MOSAIC_HPP__
#define __MOSAIC_HPP__

#include <vector>
#include <string>

namespace geotools {
	
	namespace raster {

		class Mosaic {
		private:
			void (*m_fileCallback)(float);
			void (*m_overallCallback)(float);
		public:
			void setFileCallback(void (*callback)(float));
			void setOverallCallback(void (*callback)(float));
			
			void mosaic(const std::vector<std::string> &files, const std::string &outfile, float distance, int tileSize, int threads = 1);
		};

	} // raster

} // geotools


#endif