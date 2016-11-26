#ifndef __POINTNORMALIZE_HPP__
#define __POINTNORMALIZE_HPP__

#include <string>
#include <list>

#include "util.hpp"

namespace geotools {

    namespace point {

        class PointNormalizeConfig {
        public:
            std::string terrainFile;
            std::string pointOutputDir;
            std::list<std::string> pointFiles;
            uint32_t threads;
            bool dropNegative;
            PointNormalizeConfig() :
                dropNegative(true) {
                
            }
        };

        class PointNormalize {
        public:
            void normalize(const PointNormalizeConfig &config, const geotools::util::Callbacks *callbacks = nullptr);
        };

    } // point

} // geotools

#endif
