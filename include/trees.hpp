#ifndef __TREETOPS_HPP__
#define __TREETOPS_HPP__

#include <string>
#include <map>

#include "util.hpp"

namespace geotools {

	namespace trees {

		class Trees;

		namespace config {


                        /**
                         * Contains configuration information for performing tree top extraction.
                         */
			class TreeTopConfig {
				friend class geotools::trees::Trees;
			public:
                                /**
                                 * Set to true to perform smoothing on the input raster.
                                 * This will force a check that the smoothing params are valild.
                                 * A smoothed filename must be provided for output, plus a sigma
                                 * and window size.
                                 */
                                bool performSmoothing;

                                /**
                                 * The size of the smoothing window >=3; an odd number.
                                 */
                                int smoothingWindow;

                                /**
                                 * The std. deviation used for generating the gaussian kernel.
                                 * 0 < n <= 1.
                                 */
                                double smoothingSigma;

                                /**
                                 * The size of the kernel used for searching for tops. 
                                 * >= 3; an odd number.
                                 */
                                int searchWindow;

                                /**
                                 * The minimum height of a pixel that will be considered for inclusion as a 
                                 * tree top.
                                 */
                                double minHeight;

                                /**
                                 * The input raster. If smoothing is to be done, an unsmoothed raster.
                                 * If smoothing is not to be done, input an already smoothed raster (recommended.)
                                 */
                                std::string inputFilename;

                                /**
                                 * If smoothing is to be done, this is the smoothed raster output file.
                                 */
                                std::string smoothedFilename;

				/**
                                 * The filename of the output table containing tree top points.
                                 * Must be an sqlite file. If the file exists and has the correct
                                 * structure, it will be cleared. If it does not exist, it is creaated.
                                 * Otherwise, a failure occurs.
                                 */
                                std::string outputFilename;


                                /**
                                 * Defines the boundaries of work to be performed. Every step of the process
                                 * will be confined, including smoothing and searching.
                                 */
                                geotools::util::Bounds bounds;

				/**
				 * A spatial reference ID for the output files.
				 */
				int srid;

				/**
				 * If true, build the index on the tops table. Can be slow.
				 */
				bool buildIndex;

				/**
				 * The cache size (B) for the sqlite database. A performance optimization.
				 */
				int tableCacheSize;

				/**
				 * The cache size (B) for rows when reading the raster.
				 */
				int rowCacheSize;

				/**
				 * Build a TreeTopConfig with defaults.
				 */
                                TreeTopConfig() :
                                        performSmoothing(false),
                                        smoothingWindow(-1),
                                        smoothingSigma(-1),
                                        searchWindow(3),
                                        minHeight(4),
					srid(0), 
					buildIndex(false),
					tableCacheSize(1024 * 1024), 
					rowCacheSize(24 * 1024 * 1024) {
                                }

				/**
				 * Set the table cache.
				 */
				void setTableCacheSize(int size) {
					tableCacheSize = size;
				}

				/**
				 * Set the row cache size.
				 */
				void setRowCacheSize(int size) {
					rowCacheSize = size;
				}

				/**
				 * Set true to build the table index.
				 */
				void setBuildIndex(bool build) {
					buildIndex = build;
				}

				/**
				 * Set the output SRID.
				 */
				void setSRID(int srid) {
					this->srid = srid;
				}

				/**
				 * Set the input raster file.
				 */
				void setInputFilename(const std::string &filename) {
					inputFilename.assign(filename);
				}

				/**
				 * Set the output database file (sqlite).
				 */
				void setOutputFilename(const std::string &filename) {
					outputFilename.assign(filename);
				}

				const std::string getOutputFilename() const {
					return outputFilename;
				}

				/**
				 * Set the size of the search window.
				 */
				void setSearchWindow(int window) {
					searchWindow = window;
				}

				/**
				 * Set the minimum pixel height for selection as a tree crown.
				 */
				void setMinHeight(double height) {
					minHeight = height;
				}

				/**
				 * Set the smoothing parameters.
				 */
                                void setSmoothing(const std::string &filename, double sigma = 0.8408964, int window = 3) {
                                        smoothedFilename.assign(filename);
                                        smoothingSigma = sigma;
                                        smoothingWindow = window;
                                        performSmoothing = true;
                                }

				/**
				 * Set the smoothed filename.
				 */
				void setSmoothedFilename(const std::string &filename) {
					smoothedFilename.assign(filename);
					performSmoothing = true;
				}

				/**
				 * Set the smoothing std. deviation.
				 */
				void setSmoothingSigma(double sigma) {
					smoothingSigma = sigma;
					performSmoothing = true;
				}

				/**
				 * Set the smoothing window size.
				 */
				void setSmoothingWindow(int window) {
					smoothingWindow = window;
					performSmoothing = true;
				}

				/**
				 * Check the validity of the configuration.
				 */
                                void check() const {
					if(performSmoothing) {
						if(smoothedFilename.empty())
							g_argerr("Smoothed raster filename must not be empty.");
						if(smoothingSigma <= 0 || smoothingSigma > 1)
							g_argerr("Std. deviation for smoothing must be 0 < n <= 1. " << smoothingSigma << " given.");
						if(smoothingWindow < 3)
							g_argerr("Smoothing window must be an odd number greater than or equal to 2. " << smoothingWindow << " given.");
						if(smoothingWindow % 2 == 0) 
							g_argerr("The smooothing window must be odd and >=3.");
					}
					if(inputFilename.empty())
						g_argerr("The input raster filename must not be empty.");
					if(outputFilename.empty())
						g_argerr("The output database filename must not be empty.");
					if(inputFilename == outputFilename)
						g_argerr("The input and output filenames must be different.");
					if(performSmoothing && inputFilename == smoothedFilename)
						g_argerr("The input and smoothed filenames must be different.");
					if(searchWindow < 3)
						g_argerr("Search window must be an odd number greater than or equal to 2. " << searchWindow << " given.");
					if(searchWindow % 2 == 0) 
						g_argerr("The search window must be odd and >=3.");
				}

                        };

                } // config


		namespace util {

			/**
			 * A simple class for maintaining information about a tree top.
			 */
			class Top {
			public:
				size_t id;
				double x;
				double y;
				double z;
				int col;
				int row;
				
				Top(size_t id, double x, double y, double z, int col, int row);
				
				Top(const Top &t);
				
				Top();

			};

		} // util


		class DLL_EXPORT Trees {
		private:
			/**
			 * A convenience method for smoothing the input raster before using it to generate crowns
			 * or treetops.
			 *
			 * inraster   - The raster to be smoothed.
			 * outraster  - The smoothed raster.
			 * sigma      - The standard deviation of the gaussian kernel. Default 0.8408964.
			 * window     - The size of the window to use for smoothing. Default 3.
			 */
			void smooth(const std::string &inraster, const std::string &outraster, double sigma, double window);

		public:

			/**
			 * Locates tree top points on a canopy height model.
			 *
			 * config     - A TreeTopConfig opject containing running parameters.
			 * tops       - If not null, this vector is populated with objects representing the tree tops.
			 */
			void treetops(const geotools::trees::config::TreeTopConfig &config, std::vector<std::unique_ptr<geotools::trees::util::Top> > *tops = nullptr);

			/**
			 * Performs tree crown delineation using a (preferrably smoothed) input raster and a
			 * vector file (sqlite) containing tree tops as seeds. Output is an integer raster with 
			 * cell values* representing tree top IDs, and an optional vector which is the polygonized 
			 * version* of the raster. The table should have been generated using the treetops() method
			 * to ensure that its structure is correct.
			 *
			 * inraster   - The input raster. Preferrably smoothed (i.e. using smooth()).
			 * tops       - A vector containing tops, generated by treetops().
			 * crownsrast - The raster to write the crowns to (required.)
			 * crownsvect - The vector to write the crowns to (sqlite; required.)
			 * threshold  - A ratio (0 < n < 1) which filters out pixels that are less than a given
			 *		proportion of the tree top's height.
			 * radius     - The maximum radius of an idividual tree crown.
			 * minHeight  - Heights below this value will not be considered.
			 * d8         - Use D8 search rather than D4.
			 */
			void treecrowns(const std::string &inraster, const std::string &topsvect, const std::string &crownsrast, 
				const std::string &crownsvect, double threshold, double radius, double minHeight, bool d8 = false);

		};

	} // trees

} // geotools

#endif
