/*
 * treecrowns.cpp
 *
 *  Created on: May 3, 2016
 *      Author: rob
 */

#include <queue>
#include <iostream>

#include "Raster.hpp"
#include "Vector.hpp"

static unsigned int _crown_id = 0;

/**
 * A rule provides parameters for processing crowns given
 * their properties, such as height.
 */
class Rule {
public:
	double minHeight;
	double maxHeight;
	double dropOff;
	double radius;
	Rule(double minHeight, double maxHeight, double dropOff, double radius) {
		this->minHeight = minHeight;
		this->maxHeight = maxHeight;
		this->dropOff = dropOff;
		this->radius = radius;
	}
	bool match(double height) {
		return height >= this->minHeight && height < this->maxHeight;
	}
};

/**
 * Represents a single cell (pixel) in the input/output rasters.
 */
class Px {
public:
	double crownHeight;
	int col;
	int row;
	unsigned int id;
	unsigned int radius;
	Px(double crownHeight, int col, int row, unsigned int id, unsigned int radius = 0) {
		this->crownHeight = crownHeight;
		this->col = col;
		this->row = row;
		this->id = id;
		this->radius = radius;
	}
};

/** 
 * Represents a single crown.
 */
class Crown {
public:
	unsigned int id;
	int col;
	int row;
	double height;
	double radius;
	Crown(int col, int row, double height, double radius) {
		this->col = col;
		this->row = row;
		this->height = height;
		this->radius = radius;
		this->id = ++_crown_id;
	}
};

/**
 * Returns a rule corresponding to the given crown height.
 */
// TODO: Don't use a pointer; use a reference to a special no-rule Rule.
Rule *getRuleForHeight(std::vector<Rule> &rules, double height) {
	for(Rule &rule:rules) {
		if(rule.match(height))
			return &rule;
	}
	return nullptr;
}

/**
 * Returns true if the given cell is a maximum, according to the given Rule.
 */
bool isMaximum(Raster<float> &raster, int col, int row, std::vector<Rule> &rules, double *height, double resolution) {

	*height = raster.get(col, row);
	if(*height == raster.nodata())
		return false;

	Rule *rule = getRuleForHeight(rules, *height);
	if(rule == nullptr)
		return false;

	double radius = rule->radius;
	int ks = _max(3, (int) (radius / resolution / 2.0) * 2 + 1);
	int ks0 = (int) ks / 2;

	if(col < ks0 || row < ks0 || col >= raster.cols() - ks0 || row >= raster.rows() - ks0)
		return false;

	for(int r = 0; r < ks; ++r) {
		for(int c = 0; c < ks; ++c) {
			if(r == ks0 && c == ks0) continue;
			if(raster.get(col - ks0 + c, row - ks0 + r) > *height)
				return false;
		}
	}
	return true;
}

/**
 * Find the maxima in the raster and append a Crown object to the list for each one.
 */
void findMaxima(Raster<float> &raster, std::vector<Crown> &crowns, std::vector<Rule> &rules, double resolution) {
	double height = nan("");
	for(int r = 0; r < raster.rows(); ++r) {
		for(int c = 0; c < raster.cols(); ++c) {
			if(isMaximum(raster, c, r, rules, &height, resolution))
				crowns.push_back(Crown(c, r, height, 0.0));
		}
	}
}

/**
 * Compute the table of Gaussian weights given the size of the table
 * and the std. deviation.
 */
void gaussianWeights(double *weights, int size, double sigma) {
	if(size % 2 == 0) ++size;
	for(int r = 0; r < size; ++r) {
		for(int c = 0; c < size; ++c) {
			int x = c - size / 2;
			int y = r - size / 2;
			weights[r * size + c] = (1 / (2 * PI * sigma * sigma)) * pow(E, -((x*x + y*y) / (2.0*sigma*sigma)));
		}
	}
}

/**
 * Smooth the raster and write the smoothed version to the output raster.
 */
void smooth(Raster<float> &raster, Raster<float> &smoothed, double sigma, int size) {
	double weights[size*size];
	gaussianWeights(weights, size, sigma);
	for(int r = 0; r < raster.rows(); ++r) {
		for(int c = 0; c < raster.cols(); ++c) {
			double t = 0.0;
			for(int gr = 0; gr < size; ++gr) {
				for(int gc = 0; gc < size; ++gc) {
					int rc = c - size / 2 + gc;
					int rr = r - size / 2 + gr;
					if(rc >= 0 && rc < raster.cols() && rr >= 0 && rr < raster.rows())
						t += weights[gr * size + gc] * raster.get(rc, rr);
				}
			}
			smoothed.set(c, r, t);
		}
	}
}

/**
 * Delineate and output the treetops present in the given canopy height model.
 * Outputs a tree tops vector (with heights, IDs), a smoothed CHM and a raster containing 
 * the crowns (with IDs).
 */
// TODO: Output crown vector.
void treetops(std::string &inraster, std::string &smraster, std::string &crownshp, std::string &topshp) {

	if(inraster.empty())
		throw "Input raster cannot be empty.";
	if(crownshp.empty())
		throw "The crown output name cannot be empty.";
	if(topshp.empty())
		throw "The treetop output name cannot be empty.";
	if(smraster.empty())
		throw "The smoothed raster output name cannot be empty.";
	if(smraster == topshp || topshp == crownshp || smraster == crownshp)
		throw "The output file names are in conflict.";

	std::vector<Rule> rules = {
		Rule(3.0, 6.0, 0.65, 1.0),
		Rule(6.0, 35.0, 0.65, 1.0),
		Rule(35.0, 100.0, 0.65, 1.0)
	};

	std::cerr << "Loading " << inraster << std::endl;

	// Work grid.
	static Grid<float> k(3, 3);

	// Input raster.
	Raster<float> raster(inraster);

	// Output raster.
	std::string proj;
	raster.projection(proj);

	// To store the IDs and extents of the crowns.
	Raster<int> ids(crownshp, raster.minx(), raster.miny(), raster.maxx(), raster.maxy(),
				raster.resolution(), 0, proj);
	// So track which cells are visited (faster than using IDs).
	Grid<char> visited(raster.cols(), raster.rows());
	visited.fill(0);

	// Gaussian-smoothed version of source raster.
	Raster<float> smoothed(smraster, raster);

	// Treetops vector file.
	remove(topshp.c_str());
	std::map<std::string, int> attribs;
	attribs["id"] = Vector::INTEGER;
	attribs["height"] = Vector::DOUBLE;
	Vector tops(topshp, Vector::POINT, proj, attribs);

	// Tree crown objects.
	std::vector<Crown> crowns;

	// Smooth the input raster.
	smooth(raster, smoothed, 0.84089642, 7);

	// Find the treetops.
	findMaxima(smoothed, crowns, rules, raster.resolution());

	std::cerr << "Building initial queue." << std::endl;

	// Create a queue of cells contraining tree tops (maxima).
	std::queue<Px> cells;
	for(Crown &crn:crowns) {
		cells.push(Px(crn.height, crn.col, crn.row, crn.id));
		std::unique_ptr<Geom> geom = tops.addPoint(
			raster.toX(crn.col) + raster.resolutionX() / 2.0,
			raster.toY(crn.row) + raster.resolutionY() / 2.0
		);
		geom->setAttribute("id", (int) crn.id);
		geom->setAttribute("height", crn.height);
	}

	std::cerr << "Finding crowns. " << crowns.size() << std::endl;

	// Iterate over the cell queue.
	while(!cells.empty()) {
		Px px = cells.front();
		cells.pop();
		// Ignore edge cells. Skip visited cells.
		if(px.col < 1 || px.row < 1
				|| px.col >= smoothed.cols() - 1 || px.row >= smoothed.rows() - 1
				|| visited.get(px.col, px.row))
			continue;
		// Update the output with the id of the crown.
		ids.set(px.col, px.row, px.id);
		visited.set(px.col, px.row, 1);
		// Get the rule for the current crown height (may be null).
		Rule *rule = getRuleForHeight(rules, px.crownHeight);
		if(rule == nullptr)
			continue;
		// Load the local grid.
		smoothed.readBlock(px.col - 1, px.row - 1, 3, 3, k);
		// Get the height at the center of the kernel.
		double height = k.get(1, 1);
		// Iterate over kernel.
		for(int r = 0; r < 3; ++r) {
			for(int c = 0; c < 3; ++c) {
				if(c == 1 && r ==1) continue; // Skip center.
				double height0 = k.get(c, r);
				// If the pixel height is high, or it's lower than the dropoff rate, ignore.
				if(height0 > height || height0 < (px.crownHeight * rule->dropOff))
					continue;
				// Else add to queue.
				cells.push(Px(px.crownHeight, px.col - 1 + c, px.row - 1 + r, px.id, px.radius + 1));
			}
		}
	}
}

void usage() {
	std::cerr << "Usage: treecrowns <options>\n"
			<< "This program delineates tree crowns from a canopy height model" << std::endl
			<< " -i -- The input raster; a LiDAR-derived canopy height model."  << std::endl
			<< " -t -- The treetop vector file. A shapefile."  << std::endl
			<< " -o -- The crown file. A raster containing cells set to the crown ID."  << std::endl
			<< " -s -- The smoothed crown file. This is an intermediate product, but"  << std::endl
			<< "       the user may wish to inspect it to tweak the parameters." << std::endl;
}

int main(int argc, char **argv) {

	std::string inraster; // Input raster.
	std::string crownshp; // Crowns output.
	std::string topshp;   // Treetops output.
	std::string smraster; // Smoothed raster (optional).

	int i = 1;
	for(; i < argc; ++i) {
		std::string arg = argv[i];
		if(arg == "-i") {
			inraster = argv[++i];
		} else if(arg == "-o") {
			crownshp = argv[++i];
		} else if(arg == "-t") {
			topshp = argv[++i];
		} else if(arg == "-s") {
			smraster = argv[++i];
		}
	}

	try {
		treetops(inraster, smraster, crownshp, topshp);
	} catch(const char *e) {
		std::cerr << e << std::endl;
		usage();
		return 1;
	}

	return 0;
}



