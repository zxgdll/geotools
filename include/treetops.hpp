#ifndef __TREETOPS_HPP__
#define __TREETOPS_HPP__

#include <string>

#include "geotools.h"
#include "util.hpp"

namespace geotools {

    namespace treetops {

        namespace config {

            // Contains configuration information for performing tree top extraction.

            class TreetopsConfig {
            public:

                // Defines the boundaries of work to be performed. Every step of the process
                // will be confined, including smoothing and searching.
                geotools::util::Bounds bounds;

                // A spatial reference ID for the output files.
                int srid;

                // If true, build the index on the tops table. Can be slow.
                bool buildIndex;

                // The cache size (B) for the sqlite database. A performance optimization.
                int tableCacheSize;

                // The cache size (B) for rows when reading the raster.
                int rowCacheSize;

                // The number of threads to use in execution.
                uint8_t threads;

                // Set to true to perform smoothing on the input raster.
                // This will force a check that the smoothing params are valild.
                // A smoothed filename must be provided for output, plus a sigma and window size.
                bool doSmoothing;

                // The size of the smoothing window >=3; an odd number.
                int smoothWindowSize;

                // The std. deviation used for generating the gaussian kernel.
                // 0 < n <= 1.
                double smoothStdDev;

                std::string smoothOriginalCHM;
                std::string smoothSmoothedCHM;

                // If true, treetop location will be performed.
                bool doTops;

                // The minimum height of a pixel that will be considered for inclusion as a 
                // tree top.
                double topsMinHeight;

                // The size of the top finding window >=3; an odd number.
                int topsWindowSize;

                std::string topsOriginalCHM;
                std::string topsSmoothedCHM;
                std::string topsTreetopsDatabase;

                // Set to true to delineate crowns.
                bool doCrowns;

                // The maximum crown radius.
                double crownsRadius;

                // The maximum height of a crown as a fraction of top height.
                double crownsHeightFraction;

                // The maximum height of pixels to consider for inclusion.
                double crownsMinHeight;

                // The input raster -- ideally the same one used for tops.
                std::string crownsSmoothedCHM;

                // The treetops database file.
                std::string crownsTreetopsDatabase;

                std::string crownsCrownsRaster;
                std::string crownsCrownsDatabase;

                // Build a TreetopsConfig with defaults.

                TreetopsConfig() :
                srid(0),
                buildIndex(false),
                tableCacheSize(1024 * 1024),
                rowCacheSize(24 * 1024 * 1024),
                doSmoothing(false),
                smoothWindowSize(3),
                smoothStdDev(0.8),
                doTops(false),
                topsMinHeight(4.0),
                topsWindowSize(7),
                doCrowns(false),
                crownsRadius(10.0),
                crownsHeightFraction(0.65),
                crownsMinHeight(4.0) {
                }

                void checkSmoothing() const {
                    if (!doSmoothing)
                        g_argerr("Not configured to perform smoothing.");
                    if (smoothOriginalCHM.empty())
                        g_argerr("Smoothing: CHM filename must not be empty.");
                    if (smoothSmoothedCHM.empty())
                        g_argerr("Smoothing: Output filename must not be empty.");
                    if (smoothStdDev <= 0 || smoothStdDev > 1)
                        g_argerr("Smoothing: Std. deviation must be 0 < n <= 1. " << smoothStdDev << " given.");
                    if (smoothWindowSize % 2 == 0 || smoothWindowSize < 3)
                        g_argerr("Smoothing: The window must be odd and >=3.");
                }

                void checkTops() const {
                    if (!doTops)
                        g_argerr("Not configured to find treetops.");
                    if (topsOriginalCHM.empty())
                        g_argerr("Tops: Unsmoothed CHM filename must not be empty.");
                    if (topsSmoothedCHM.empty())
                        g_argerr("Tops: Smoothed CHM filename must not be empty.");
                    if (topsTreetopsDatabase.empty())
                        g_argerr("Tops: Treetops database filename must not be empty.");
                    if (topsWindowSize % 2 == 0 || topsWindowSize < 3)
                        g_argerr("Tops: Treetops window size must be an odd number >= 3. " << topsWindowSize << " given.");
                }

                void checkCrowns() const {
                    if (!doCrowns)
                        g_argerr("Not configured to find crowns.");
                    if (crownsRadius <= 0.0)
                        g_argerr("Crowns: The maximum crown radius must be > 0. " << crownsRadius << " given.");
                    if (crownsHeightFraction <= 0.0 || crownsHeightFraction > 1.0)
                        g_argerr("Crowns: The crown height fraction must be between 0 and 1. " << crownsHeightFraction << " given.");
                    if (crownsCrownsRaster.empty())
                        g_argerr("Crowns: Output raster filename must not be empty.")
                        if (crownsTreetopsDatabase.empty())
                            g_argerr("Crowns: Treetops database filename must not be empty.");
                    if (crownsSmoothedCHM.empty())
                        g_argerr("Crowns: Smoothed CHM filename must not be empty.");
                }

                // Check the validity of the configuration.

                void check() const {
                    if (doSmoothing)
                        checkSmoothing();
                    if (doTops)
                        checkTops();
                    if (doCrowns)
                        checkCrowns();
                }

            };

        } // config

        namespace util {

            // A simple class for maintaining information about a tree top.

            class Top {
            public:
                uint64_t id;
                double x;
                double y;
                double z;
                double uz; // The unsmoothed z value.
                int col;
                int row;

                Top(uint64_t id, double x, double y, double z, double uz, int col, int row);

                Top();
            };

        } // util

        class DLL_EXPORT Treetops {
        private:
            geotools::util::Callbacks *m_callbacks;

        public:

            void setCallbacks(geotools::util::Callbacks *callbacks);

            // A convenience method for smoothing the input raster before using it to generate crowns
            // or treetops.
            void smooth(const geotools::treetops::config::TreetopsConfig &config);

            // Locates tree top points on a canopy height model.
            void treetops(const geotools::treetops::config::TreetopsConfig &config);

            // Performs tree crown delineation using a (preferrably smoothed) input raster and a
            // vector file (sqlite) containing tree tops as seeds. Output is an integer raster with 
            // cell values representing tree top IDs, and an optional vector which is the polygonized 
            // version of the raster. The table should have been generated using the treetops() method
            // to ensure that its structure is correct.
            void treecrowns(const geotools::treetops::config::TreetopsConfig &config);

        };

    } // trees

} // geotools

#endif
