#ifndef __SPECTRAL_HPP__
#define __SPECTRAL_HPP__

#include <vector>
#include <string>
#include <set>

namespace geotools {

        namespace spectral {

		namespace config {

			class SpectralConfig {
			public:
				std::vector<std::string> spectralFilenames;
				std::string indexFilename;
				std::string outputFilename;
				std::set<int> bands;
				unsigned short nodata;
			};


		} // config

                class Spectral {
		public:
			void extractSpectra(const geotools::spectral::config::SpectralConfig &config);

                };

        } // spectral

} // geotools


#endif                                                                                                                               
