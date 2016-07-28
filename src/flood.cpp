/*****************************************************************************************************
 * This program 'floods' a digital terrain model using the flood fill algorithm, generating raster
 * and vector maps of the basin extents and 'spill points' indicating locations where basins can 
 * be expected to connect. The user can provide seed points to the program, or the program
 * can generate random seeds. The user can provide start and stop elevations, otherwise the program
 * will start at the elevation of the lowest seed, and stop at the elevation above the highest
 * seed when all basins have converged to one.
 *
 * The output of this program can have many uses, such as informing the selection of outlets
 * for watershed mapping with r.water.outlet.
 *s
 * Run the program with no arguments for instructions.
 *****************************************************************************************************/

#include <vector>
#include <ostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <sstream>

#include "Raster.hpp"

namespace flood {

	// Reads a line from a CSV stream into a vector of values.
	bool readline(std::istream &st, std::vector<std::string> &row) {
		row.clear();
	    std::string line;
	    std::getline(st, line);
	    std::stringstream linest(line);
	    std::string cell;
	    bool data = false;
	    while(std::getline(linest, cell, ',')) {
	        row.push_back(cell);
	        data = true;
	    }
	    return data;
	}

	// Represents a bounding box with columns and rows.
	class Box {
	private:
		int m_c0, m_c1, m_r0, m_r1;
	public:
		Box(int c0, int r0, int c1, int r1) :
			m_c0(c0), m_r0(r0),
			m_c1(c1), m_r1(r1) {
		}
		Box() : Box(0, 0, 0, 0) {
		}
	};

	static unsigned int __cell_id = 0;

	// Represents a grid cell.
	class Cell {
	private:
		int m_col;
		int m_row;
		unsigned int m_id;
		float m_value;

	public:
		Cell(int col, int row) : 
			m_col(col),
			m_row(row) {
			m_id = ++__cell_id;
		}

		Cell(unsigned int id, int col, int row) :
			m_id(id), 
			m_col(col), m_row(row) {
			__cell_id = _max(id, __cell_id);
		}

		Cell(unsigned int id, int col, int row, float value) :
			m_id(id), 
			m_col(col), m_row(row), 
			m_value(value) {
			__cell_id = _max(id, __cell_id);

		}

		Cell(int col, int row, float value) :
			m_col(col), m_row(row), 
			m_value(value) {
			m_id = ++__cell_id;

		}

		int row() const {
			return m_row;
		}

		int col() const {
			return m_col;
		}

		unsigned int id() const {
			return m_id;
		}

		float value() const {
			return m_value;
		}

		double distance(Cell &other, double resx, double resy) const {
			double x0 = col() * resx;
			double y0 = row() * resy;
			double x1 = other.col() * resx;
			double y1 = other.row() * resy;
			return std::sqrt(_sq(x0 - x1) + _sq(y0 - y1));
		}

	};

	// Represents a basin, with bounds and area and the ID of the seed that spawned it.
	class Basin {
	private:
		int m_minc;
		int m_minr;
		int m_maxc;
		int m_maxr;
		int m_area;
		unsigned int m_id;

	public:
		Basin(unsigned int id, int minc, int minr, int maxc, int maxr, int area) :
			m_id(id),
			m_minc(minc), m_minr(minr), 
			m_maxc(maxc), m_maxr(maxr),
			m_area(area) {

		}

		// Make a list of all the cells that are on the edges of this basin.
		std::vector<Cell> computeEdges(Grid<unsigned int> &grd) {
			std::vector<Cell> edgeCells;
			for(int r = m_minr; r < m_maxr; ++r) {
				for(int c = m_minc; c < m_maxc; ++c) {
					bool edge = false;
					for(int rr = r - 1; !edge && rr < r + 2; ++rr) {
						for(int cc = c - 1; !edge && cc < c + 2; ++cc) {
							if(grd.has(cc, rr) && grd.get(cc, rr) != m_id) {
								edgeCells.push_back(Cell(c, r));
								edge = true;
							}
						}
					}
				}
			}
			return edgeCells;
		}

	};

	// A flood fill operator that fills pixels whose values are lower than
	// the given elevation.
	template <class T>
	class LEFillOperator : public FillOperator<T> {
	private:
		T m_elevation;
	public:
		LEFillOperator(T elevation) : 
			m_elevation(elevation) {
		}
		bool fill(T value) const {
			return value <= m_elevation;
		}
	};

	// Represents a spill point between two cells.
	class SpillPoint {
	private:
		Cell m_c1;
		Cell m_c2;
	public:
		SpillPoint(Cell &c1, Cell &c2) :
			m_c1(c1), m_c2(c2) {
		}
		const Cell& cell1() const {
			return m_c1;
		}
		const Cell& cell2() const {
			return m_c2;
		}
	};

	class Config {
	public:
		int m_t; // number of threads
		double m_start;
		double m_end;
		double m_step;
		double m_minBasinArea;
		double m_maxSpillDist;
		std::string m_vdir;
		std::string m_rdir;
		std::string m_spill;
		std::string m_input;
		Raster<float> m_dem;
		Raster<unsigned int> m_basins;
		std::vector<Cell> m_seeds;
		std::vector<Basin> m_basinList;
		std::vector<SpillPoint> m_spillPoints;

		Config(std::string &input, std::string &vdir, std::string &rdir, std::string &spill,
			double start, double end, double step, double minBasinArea, double maxSpillDist) :
			m_input(input), 
			m_vdir(vdir), m_rdir(rdir), m_spill(spill),
			m_start(start), m_end(end), m_step(step),
			m_minBasinArea(minBasinArea), m_maxSpillDist(maxSpillDist) {
		}

		~Config() {
		}

		double minBasinArea() {
			return m_minBasinArea;
		}

		double maxSpillDist() {
			return m_maxSpillDist;
		}

		void loadSeeds(std::string &seedFile, bool header = false) {
			if(seedFile.empty())
				_argerr("No seed file given.");
			std::ifstream csv(seedFile);
			std::vector<std::string> row;
			if(header) 
				readline(csv, row);
			while(readline(csv, row)) {
				unsigned int id = atoi(row[0].c_str());
				double x = atof(row[1].c_str());
				double y = atof(row[2].c_str());
				m_seeds.push_back(Cell(id, m_dem.toCol(x), m_dem.toRow(y)));
			}
			_trace("Seeds loaded from file: " << m_seeds.size());
		}

		std::vector<Cell>& seeds() {
			return m_seeds;
		}

		Grid<float>& dem() {
			return m_dem;
		}

		Grid<unsigned int>& basins() {
			return m_basins;
		}

		void init() {
			_trace("Checking...");
			if(m_input.empty())
				_argerr("Input DEM must be provided.");
			if(m_vdir.empty() && m_rdir.empty() && m_spill.empty())
				_argerr("At least one of vector directory, raster directory, or spill point file must be given.");
			if(m_vdir.empty())
				_trace("WARNING: No vector directory; not producing flood polygons.");
			if(m_rdir.empty())
				_trace("WARNING: No raster directory; not producing flood rasters.");
			if(m_spill.empty())
				_trace("WARNING: No raster directory; not producing spill points.");
			if(std::isnan(m_start))
				_trace("WARNING: Start value not given; using raster minimum.");
			if(std::isnan(m_end))
				_trace("WARNING: End value not given; using raster maximum.");
			if(m_step <= 0.0)
				_argerr("The step elevation must be greater than zero.");
			if(m_t < 1)
				_trace("WARNING: Invalid number of threads. Using 1.");
			if(m_minBasinArea <= 0.0)
				_argerr("Min basin area must be greater than zero.");
			if(m_maxSpillDist <= 0.0)
				_argerr("Max spill distance must be greater than zero.");

			_trace("Initing...");
			m_dem.init(m_input);
			m_basins.init("basins.tif", m_dem);
			m_basins.nodata(0);

			if(std::isnan(m_start)) {
				m_start = m_dem.min();
			} else {
				m_start = _max(m_dem.min(), m_start);
			}
			if(std::isnan(m_end)) {
				m_end = m_dem.max();
			} else {
				m_end = _min(m_dem.max(), m_end);
			}
			if(m_end <= m_start)
				_argerr("The ending elevation must be larger than the starting elevation.");
		}

		/**
		 * Perform flood filling and identify basins.
		 */
		int fillBasins(float elevation) {
			_trace("Filling basins.");

			std::string filename;
			if(!m_rdir.empty()) {
				std::stringstream fn;
				fn << m_rdir << "/" << "basin_" << std::setfill('0') 
					<< std::setw(3) << elevation << ".tif";
				filename.assign(fn.str());
			} else {
				filename.assign("/tmp/basin.tif");
			}

			m_basins.fill(0);
			m_basinList.clear();
			
			LEFillOperator<float> op(elevation);

			for(Cell seed : seeds()) {

				_trace("Seed: " << seed.id() << ": " << seed.col() << ", " << seed.row());

				if(!m_dem.has(seed.col(), seed.row())) {
					_trace("WARNING: Found a seed out of bounds.");
					continue;
				}

				// Fill the basin based on the elevations in dem.
				std::vector<int> result = m_dem.floodFill(seed.col(), seed.row(), op, m_basins, seed.id());
				int area = result[4];

				_trace("Basin: area: " << area);

				if(area >= minBasinArea()) {
					// If it's large enough, save the basin.
					m_basinList.push_back(Basin(seed.id(), result[0], result[1], result[2], result[3], area));
				} else {
					// If the basin is too small, fill it with nodata. Do not collect more spill points.
					m_basins.floodFill(seed.col(), seed.row(), seed.id(), m_basins.nodata());
				}

			}
			
			return m_basinList.size();
		}

		void saveBasinRaster(double elevation) {
			std::stringstream ss;
			ss << m_rdir << "/" << (int) (elevation / m_step) << ".tif";
			_trace("Saving: " << ss.str());
			Raster<unsigned int> r(ss.str(), m_basins);
			r.writeBlock(m_basins);
		}

		void saveBasinVector() {

		}

		bool findSpillPoints() {
			_trace("Finding spill points.");

			m_spillPoints.clear();

			// Compare each basin to each other basin.
			for(int i = 0; i < m_basinList.size(); ++i) {
				for(int j = i + 1; j < m_basinList.size(); ++j) {

					Basin b0 = m_basinList[i];
					Basin b1 = m_basinList[j];
					std::vector<flood::Cell> cells0 = b0.computeEdges(m_basins);
					std::vector<flood::Cell> cells1 = b1.computeEdges(m_basins);

					// Compare the distances; save the ones that are near enough.
					for(int k = 0; k < cells0.size(); ++k) {
						for(int l = 0; l < cells1.size(); ++l) {
							if(cells0[k].distance(cells1[l], m_dem.resolutionX(), m_dem.resolutionY()) <= m_maxSpillDist)
								m_spillPoints.push_back(SpillPoint(cells0[k], cells1[l]));
						}
					}
				}
			}

			return m_spillPoints.size();
		}

		// Output the spill points to a stream, with comma delimiters.
		// The fields are: ID1, x1, y1, ID2, x2, y2, midpoint x, midpoint y, distance
		void saveSpillPoints(std::ostream &out) {
			_trace("Outputting spill points.");

			for(SpillPoint sp : m_spillPoints) {
				const flood::Cell &c1 = sp.cell1();
				const flood::Cell &c2 = sp.cell2();
				double x1 = c1.col() * m_dem.resolutionX() + m_dem.minx();
				double y1 = c1.row() * m_dem.resolutionY() + m_dem.maxy();
				double x2 = c2.col() * m_dem.resolutionX() + m_dem.minx();
				double y2 = c2.row() * m_dem.resolutionY() + m_dem.maxy();
				double x3 = (x1 + x2) / 2.0;
				double y3 = (y1 + y2) / 2.0;
				double dist = std::sqrt(_sq(x1 - x2) + _sq(y1 - y2));
				out << c1.id() << "," << x1 << "," << y1 << "," << c2.id() << "," << x2 << "," << y2 << "," << x2 << "," << y3 << "," << dist << std::endl;
			}

		}

		/**
		 * Find the cells at the bottoms of depressions.
		 */
		void findMinima() {
			_trace("Finding minima.");

			m_seeds.clear();
			for(int r = 0; r < m_dem.rows(); ++r) {
				for(int c = 0; c < m_dem.cols(); ++c) {
					bool skip = false;
					if(m_dem.isNoData(c, r)) 
						continue;
					for(int rr = _max(0, r - 1); !skip && rr < _min(r + 2, m_dem.rows()); ++rr) {
						for(int cc = _max(0, c - 1); !skip && cc < _min(c + 2, m_dem.cols()); ++cc) {
							if((cc == c && rr == r) || m_dem.isNoData(cc, rr)) 
								continue;
							if(m_dem.get(cc, rr) < m_dem.get(c, r))
								skip = true;
						}
					}
					if(!skip)
						m_seeds.push_back(flood::Cell(c, r, m_dem.get(c, r)));
				}
			}
			_trace("Seeds found from minima: " << m_seeds.size());
		}
	};


	/**
	 */
	void flood(std::string &input, std::string &seeds, std::string &vdir, 
		std::string &rdir, std::string &spill, double start, double end, double step, 
		double minBasinArea, double maxSpillDist, int t) {

		_trace("Flooding...");

		// Build the config object.
		flood::Config config(input, vdir, rdir, spill, start, end, step, minBasinArea, maxSpillDist);

		config.init();

		_trace("Building seed list...");
		if(!seeds.empty()) {
			// Load the seeds if given.
			config.loadSeeds(seeds);
		} else {
			// Otherwise find and use minima.
			config.findMinima();
		}

		double elevation = start;
		while(elevation <= end) {
			_trace("Filling to " << elevation);
			if(config.fillBasins(elevation) > 0) {
				if(!rdir.empty())
					config.saveBasinRaster(elevation);
				if(!vdir.empty())
					config.saveBasinVector();
				if(config.findSpillPoints() > 0)
					config.saveSpillPoints(std::cout);
			}
			elevation += step;
		}

	}
}



void usage() {
	std::cerr << "Usage: flood <options>\n"
		<< " -i     Input elevation raster.\n"
		<< " -s     Shapefile or CSV containing seed points. If not given, minima are used.\n"
		<< " -v     Directory for basin vectors. If not given, they are not produced.\n"
		<< " -r     Directory for basin rasters. If not given, they are produced.\n"
		<< " -p     Spill point file. If not given, they are not produced. A Shapefile or CSV.\n"
		<< " -start Starting elevation, in same units as input raster.\n"
		<< " -end   Ending elevation.\n"
		<< " -step  Step elevation.\n"
		<< " -t     Number of threads to use. Default 1.\n"
		<< " -b     Minimum basin area.\n"
		<< " -d     Maximum spill distance.\n";
}

int main(int argc, char **argv) {

	_loglevel(LOG_TRACE);

	std::string input;
	std::string seeds;
	std::string vdir;
	std::string rdir;
	std::string spill;
	double start = 0.0;
	double end = 0.0;
	double step = 0.0;
	double maxSpillDist = 1.0;
	double minBasinArea = 100.0;
	int t = 1;

	for(int i = 1; i < argc; ++i) {
		std::string a(argv[i]);
		if(a == "-i") {
			input.assign(argv[++i]);
		} else if(a == "-s") {
			seeds.assign(argv[++i]);
		} else if(a == "-v") {
			vdir.assign(argv[++i]);
		} else if(a == "-r") {
			rdir.assign(argv[++i]);
		} else if(a == "-p") {
			spill.assign(argv[++i]);
		} else if(a == "-start") {
			start = atof(argv[++i]);
		} else if(a == "-end") {
			end = atof(argv[++i]);
		} else if(a == "-step") {
			step = atof(argv[++i]);
		} else if(a == "-t") {
			t = atoi(argv[++i]);
		} else if(a == "-b") {
			minBasinArea = atof(argv[++i]);
		} else if(a == "-d") {
			maxSpillDist = atof(argv[++i]);
		}
	}

	try {

		flood::flood(input, seeds, vdir, rdir, spill, start, end, step, minBasinArea, maxSpillDist, t);

	} catch(const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;
}
	
