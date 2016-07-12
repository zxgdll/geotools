/*
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

/**
 * Represents a single polygon in the input shape file.
 * An instance of this class can comute all the necessary
 * statistics for its polygon.
 */
class Stat {

private:
	// The list of values.
	std::vector<double> m_values;
	// Class.
	int m_class;
	bool m_sorted;
public:

	/**
	 * Construct a Stat object for the point class
	 * using the given geometry and feature.
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

	double min() {
		if(count() > 0) {
			sort();
			return m_values[0];
		}
		return nan("");
	}

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


void rasterstats(std::string &clsfile, std::vector<std::string> &files) {

	if(clsfile.empty())
		throw std::invalid_argument("A classification raster is required.");

	if(files.size() == 0)
		throw std::invalid_argument("At least one data raster is required.");

	std::map<std::string, std::map<std::string, std::map<int, Stat> > > stats;
	std::set<int> classes;

	Raster<char> clsrast(clsfile);

	for(int f0 = 0; f0 < files.size(); ++f0) {
		for(int f1 = f0 + 1; f1 < files.size(); ++f1) {

			Raster<float> frast0(files[f0]);
			Raster<float> frast1(files[f1]);

			if(frast0.resolutionX() != frast1.resolutionX() || frast0.resolutionY() != frast1.resolutionY())
				_runerr("Rasters must have the same resolutions.");

			int f0nodata = frast0.nodata();
			int f1nodata = frast1.nodata();

			double minx = _max(frast0.minx(), frast1.minx());
			double miny = _max(frast0.miny(), frast1.miny());
			double maxx = _min(frast0.maxx(), frast1.maxx());
			double maxy = _min(frast0.maxy(), frast1.maxy());

			int cols = (int) (maxx - minx) / frast0.resolutionX() + 1;
			int rows = (int) (miny - maxy) / frast0.resolutionY() + 1;

			//std::cerr << cols << " " << rows << std::endl;

			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c) {

					double x = c * frast0.resolutionX() + minx;
					double y = r * frast0.resolutionY() + miny;

					int c0 = frast0.toCol(x);
					int r0 = frast0.toRow(y);
					int c1 = frast1.toCol(x);
					int r1 = frast1.toRow(y);
					int cc = clsrast.toCol(x);
					int cr = clsrast.toRow(y);

					//std::cerr << c0 << " " << c1 << " " << cc << " " << r0 << " " << r1 << " " << cr << std::endl;
					if(c0 < 0 || c1 < 0 || cc < 0 || r0 < 0 || r1 < 0 || cr < 0 ||
						c0 >= frast0.cols() || c1 >= frast1.cols() || cc >= clsrast.cols() ||
						r0 >= frast0.rows() || r1 >= frast1.rows() || cr >= clsrast.rows())
						continue;

					int cls = clsrast.get(cc, cr);
					classes.insert(cls);
					stats[files[f0]][files[f1]][cls].add(frast0.get(c0, r0) - frast1.get(c1, r1));

				}
			}
		}

	}

	std::cout << "file";
	for(int cls : classes)
		std::cout << "," << cls << "_sum," << cls << "_count," << cls << "_mean," << cls << "_var," << cls << "_stddev";
	std::cout << std::endl;

	for(auto it = stats.begin(); it != stats.end(); ++it) {
		std::cout << it->first;
		for(auto it0 = it->second.begin(); it0 != it->second.end(); ++it0) {
			std::cout << "," << it0->first;
			for(int cls : classes) {
				if(it0->second.find(cls) == it0->second.end()) {
					std::cout << ",,,,,";
				} else {
					std::cout << "," << it0->second[cls].sum() << "," << it0->second[cls].count() << "," 
						<< it0->second[cls].mean() << "," << it0->second[cls].variance() << "," << it0->second[cls].stddev();
				}
			}
			std::cout << std::endl;
		}
	}
}

// TODO: If shapefile not given, computes stats for the entire point cloud.
void usage() {
	std::cerr << "Produces a table containing statistics for each class in the class raster\n"
		<< "for every pixel in the other rasters.\n"
		<< "Usage: lasstats [options] <class raster> <raster [raster [...]]>\n";
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
		
		rasterstats(clsfile, files);

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
