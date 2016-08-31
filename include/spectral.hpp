#ifndef __SPECTRAL_HPP__
#define __SPECTRAL_HPP__

#include <vector>
#include <string>
#include <set>

#include "geotools.h"

namespace geotools {

        namespace spectral {

			namespace config {

				class SpectralConfig {
				public:
					std::vector<std::string> spectralFilenames;
					std::string indexFilename;
					std::string outputFilename;
					std::set<int> bands;
					unsigned short specNodata;
					unsigned int idxNodata;
					bool hasIdxNodata;
					bool hasSpecNodata;

					SpectralConfig() :
						specNodata(0),
						idxNodata(0) {
					}

					void check() const {
						if(spectralFilenames.size() == 0)
							g_argerr("There must be at least one spectral file.");
						if(indexFilename.empty())
							g_argerr("The index filename must be given.");
						if(outputFilename.empty())
							g_argerr("The output filename must be given.");
						if(bands.size() == 0)
							g_warn("No bands given. Using all bands.");
					}

				};


			} // config

            class Spectral {
			public:
				void extractSpectra(const geotools::spectral::config::SpectralConfig &config);

            };

        } // spectra

} // geotools


#endif                                                                                                                               
