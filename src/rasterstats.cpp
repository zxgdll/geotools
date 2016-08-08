/*
 * This progam computes some basic statistics for differences between
 * pairs of rasters, aggregated by class.
 *
 * The classes come from a classification raster, which is the first argument.
 * Then, each of the given files is differenced with each other file
 * to produce pairwise mean, variance and standard deviation within
 * the class.
 * 
 *  Author: Rob Skelly rob@dijital.ca
 */

#include <iostream>
#include <cmath>
#include <map>
#include <set>
#include <map>

#include "geotools.h"
#include "Util.hpp"
#include "Raster.hpp"

using namespace geotools::util;
using namespace geotools::raster;

namespace geotools {

	namespace raster {

		namespace util {

			/**
			 * Represents a single class, for a single pair of files.
			 */
			class Stat {

			private:
				// The list of values.
				std::vector<double> m_values;
				// Class.
				int m_class;
				// True if the list is sorted.
				bool m_sorted;

			public:

				/**
				 * Construct a Stat object.
				 */
				Stat() :
					m_class(0),
					m_sorted(false) {
				}

				/**
				 * Sorts the values.
				 */
				void sort() {
					if(!m_sorted) {
						std::sort(m_values.begin(), m_values.end());
						m_sorted = true;
					}
				}

				/**
				 * Add a point value to the list.
				 * The list ceases to be sorted.
				 */
				void add(double value) {
					m_values.push_back(value);
					m_sorted = false;
				}

				/**
				 * Return the number of values.
				 */
				int count() {
					return m_values.size();
				}

				/**
				 * Return the minimum value.
				 */
				double min() {
					if(count() > 0) {
						sort();
						return m_values[0];
					}
					return nan("");
				}

				/**
				 * Return the maximum value.
				 */
				double max() {
					if(count() > 0) {
						sort();
						return m_values[m_values.size() - 1];
					}
					return nan("");
				}

				/**
				 * Return the sum of values.
				 */
				double sum() {
					double sum = 0;
					for(int i = 0; i < m_values.size(); ++i)
						sum += m_values[i];
					return sum;
				}

				/**
				 * Return the mean of values.
				 */
				double mean () {
					if(count() > 0)
						return sum() / count();
					return nan("");
				}

				/**
				 * Return the median of values.
				 */
				double median() {
					int num;
					if((num = count()) > 0) {
						sort();
						if(num % 2 == 0) {
							return (m_values[num / 2] + m_values[num / 2 - 1]) / 2.0;
						} else {
							return m_values[num / 2];
						}
					}
					return nan("");
				}

				/**
				 * Return the variance of values.
				 */
				double variance() {
					if(count() > 0) {
						double m = mean();
						double s = 0.0;
						for(int i = 0; i < m_values.size(); ++i)
							s += _sq(m_values[i] - m);
						return s / (count() - 1.0);
					}
					return nan("");
				}

				/**
				 * Return the standard deviation of values.
				 */
				double stddev() {
					if(count() > 0)
						return sqrt(variance());
					return nan("");
				}

				~Stat() {
				}


			};

		} // util

		/**
		 * Compute statistics on the pairs of files in the files vector, for
		 * each class in the classification file. Print out a CSV table.
		 */
		void rasterstats(std::string &clsfile, std::vector<std::string> &files) {

			using namespace geotools::raster::util;
			
			if(clsfile.empty())
				throw std::invalid_argument("A classification raster is required.");

			if(files.size() < 2)
				throw std::invalid_argument("At least two data rasters are required.");

			std::map<std::string, std::map<std::string, std::map<int, Stat> > > stats;
			std::set<int> classes;

			Raster<unsigned char> clsrast(clsfile);

			// Loop over every pair of files.
			for(int f0 = 0; f0 < files.size(); ++f0) {
				for(int f1 = f0 + 1; f1 < files.size(); ++f1) {

					Raster<float> frast0(files[f0]);
					Raster<float> frast1(files[f1]);

					if(frast0.resolutionX() != frast1.resolutionX() || frast0.resolutionY() != frast1.resolutionY())
						_runerr("Rasters must have the same resolution.");

					int f0nodata = frast0.nodata();
					int f1nodata = frast1.nodata();

					double minx = _max(frast0.minx(), frast1.minx());
					double miny = _max(frast0.miny(), frast1.miny());
					double maxx = _min(frast0.maxx(), frast1.maxx());
					double maxy = _min(frast0.maxy(), frast1.maxy());

					int cols = (int) (maxx - minx) / frast0.resolutionX();
					int rows = (int) (miny - maxy) / frast0.resolutionY();

					//std::cerr << files[f0] << ": " << frast0.minx() << ", " << frast0.miny() << ", " << frast0.maxx() << ", " << frast0.maxy() << std::endl;
					//std::cerr << files[f1] << ": " << frast1.minx() << ", " << frast1.miny() << ", " << frast1.maxx() << ", " << frast1.maxy() << std::endl;
					//std::cerr << "range: " << minx << ", " << miny << ", " << maxx << ", " << maxy << std::endl;

					int stop = 1;

					for(double y = miny; y < maxy; y += std::abs(frast0.resolutionY())) {
						for(double x = minx; x < maxx; x += std::abs(frast0.resolutionX())) {

							int c0 = frast0.toCol(x);
							int r0 = frast0.toRow(y);
							int c1 = frast1.toCol(x);
							int r1 = frast1.toRow(y);
							int cc = clsrast.toCol(x);
							int cr = clsrast.toRow(y);

							// If any of the indices are out of bands, skip.
							if(c0 < 0 || c1 < 0 || cc < 0 || r0 < 0 || r1 < 0 || cr < 0 ||
								c0 >= frast0.cols() || c1 >= frast1.cols() || cc >= clsrast.cols() ||
								r0 >= frast0.rows() || r1 >= frast1.rows() || cr >= clsrast.rows())
								continue;

							float v0 = frast0.get(c0, r0);
							float v1 = frast1.get(c1, r1);

							// If either pixel is nodata, skip.
							if(v0 == f0nodata || v1 == f1nodata)
								continue;

							// Add the difference to the Stat object for the pair/class
							int cls = clsrast.get(cc, cr);
							classes.insert(cls);
							stats[files[f0]][files[f1]][cls].add(v0 - v1);
						}
					}
				}

			}

			// Print the header.
			std::cout << "file1,file2,class,count,sum,min,max,mean,variance,stddev" << std::endl;

			// Print the data rows.
			for(auto it = stats.begin(); it != stats.end(); ++it) {
				for(auto it0 = it->second.begin(); it0 != it->second.end(); ++it0) {
					for(int cls : classes) {
						std::cout << it->first << "," << it0->first << "," << cls;
						if(it0->second.find(cls) == it0->second.end()) {
							std::cout << ",,,,,,,";
						} else {
							std::cout << "," << it0->second[cls].sum() << "," << it0->second[cls].count() << "," 
								<< it0->second[cls].min() << "," << it0->second[cls].max() << ","
								<< it0->second[cls].mean() << "," << it0->second[cls].variance() << "," 
								<< it0->second[cls].stddev();
						}
						std::cout << std::endl;
					}
				}
			}
		}

	} // raster

} // geotools

void usage() {
	std::cerr << "Usage: lasstats [options] <class raster> <raster [raster [...]]>\n"
		<< "	Produces a table containing statistics for each class in the class raster\n"
		<< "	for the difference in every pixel for each pair of rasters.\n";
}

int main(int argc, char ** argv) {

	std::vector<std::string> files;
	std::string clsfile;

	for(int i = 1;i < argc; ++i) {
		std::string arg(argv[i]);
		if(i == 1) {
			clsfile.assign(argv[i]);
		} else {
			files.push_back(argv[i]);
		}
	}


	try {
		
		geotools::raster::rasterstats(clsfile, files);

	} catch(const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	} catch(const char *e) {
		std::cerr << e << std::endl;
		usage();
		return 1;
	}

	return 0;
}
