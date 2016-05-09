/*
 * treecrowns.cpp
 *
 *  Created on: May 3, 2016
 *      Author: rob
 */

#include <queue>
#include <iostream>

#include "Raster.hpp"

static unsigned int _crown_id = 0;

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

int crownReverseSort(Crown &a, Crown &b) {
	if(a.height > b.height) {
		return -1;
	} else if(a.height < b.height) {
		return 1;
	}
	return 0;
}

bool isMaximum(Grid<float> &raster, int col, int row, std::vector<Rule> &rules, double *height, double resolution) {

	*height = raster.get(col, row);
	if(*height == raster.nodata())
		return false;

	double dropOff = 0.0;
	double radius = 0.0;
	bool found = false;
	for(Rule &rule:rules) {
		if(rule.match(*height)) {
			dropOff = rule.dropOff;
			//radius = rule.radius;
			found = true;
			break;
		}
	}
	if(!found) return false;

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

void findMaxima(Grid<float> &raster, std::vector<Crown> &crowns, std::vector<Rule> &rules, double resolution) {
	double height = nan("");
	for(int r = 0; r < raster.rows(); ++r) {
		for(int c = 0; c < raster.cols(); ++c) {
			if(isMaximum(raster, c, r, rules, &height, resolution))
				crowns.push_back(Crown(c, r, height, 0.0));
		}
	}
}

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

void smooth(Grid<float> &raster, double sigma, int size) {
	Grid<float> smoothed(raster.cols(), raster.rows());
	smoothed.writeBlock(raster);
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
	raster.writeBlock(smoothed);
}

void treetops(std::string &inraster, std::string &crownshp, std::string &topshp, double dropoff) {

	std::vector<Rule> rules = {
		Rule(3.0, 6.0, 0.65, 5.0),
		Rule(6.0, 35.0, 0.65, 5.0),
		Rule(35.0, 100.0, 0.65, 5.0)
	};

	std::cerr << "Loading " << inraster << std::endl;
	static Grid<float> k(3, 3);
	Raster<float> src(inraster);
	Grid<float> raster(src.cols(), src.rows());
	src.readBlock(raster);
	Grid<int> ids(raster.cols(), raster.rows());
	ids.fill(0);


	std::vector<Crown> crowns;
	smooth(raster, 0.84089642, 7);
	//std::string smoo = "../data/smoothed.tif";
	//Raster<float> smoothed(smoo, src);
	//smoothed.writeBlock(raster);
	findMaxima(raster, crowns, rules, src.resolution());

	std::cerr << "Building initial queue." << std::endl;

	// A map to keep track of the number of cells linked to a given crown ID.
	// When the count falls to zero, the shape can be output, and the entry deleted.
	//std::map<unsigned int,  unsigned int> cellCounts;

	std::queue<Px> cells;
	for(Crown &crn:crowns) {
		cells.push(Px(crn.height, crn.col, crn.row, crn.id));
		//cellCounts[crn.id] = 1;
	}

	std::cerr << "Finding crowns. " << crowns.size() << std::endl;

	while(!cells.empty()) {
		Px px = cells.front();
		cells.pop();
		if(px.col < 1 || px.row < 1 || px.col >= raster.cols() - 1 || px.row >= raster.rows() - 1)
			continue;
		if(ids.get(px.col, px.row) != 0)
			continue;
		raster.readBlock(px.col - 1, px.row - 1, 3, 3, k);
		double height = k.get(1, 1);
		ids.set(px.col, px.row, px.id);
		for(int r = 0; r < 3; ++r) {
			for(int c = 0; c < 3; ++c) {
				if(c == 1 && r ==1) 
					continue;
				double height0 = k.get(c, r);
				if(height0 > height 
					|| height0 < (px.crownHeight * dropoff)) 
					continue;
				cells.push(Px(px.crownHeight, px.col - 1 + c, px.row - 1 + r, px.id, px.radius + 1));
				//cellCounts[px.id]++;
			}
		}
		//if(cellCounts[px.id] == 0) {

		//}
	}


	std::string proj;
	src.projection(proj);
	Raster<int> output(crownshp, src.minx(), src.miny(), src.maxx(), src.maxy(),
				src.resolution(), 0, proj);
	output.writeBlock(ids);

}

int main(int argc, char **argv) {

	std::string inraster;
	std::string crownshp;
	std::string topshp;
	double dropoff = 0.65; // Percentage of elevation to drop without stop.

	int i = 1;
	for(; i < argc; ++i) {
		std::string arg = argv[i];
		if(arg == "-d") {
			dropoff = atof(argv[++i]);
		} else if(arg == "-i") {
			inraster = argv[++i];
		} else if(arg == "-o") {
			crownshp = argv[++i];
		} else if(arg == "-t") {
			topshp = argv[++i];
		}
	}

	//try {
		treetops(inraster, crownshp, topshp, dropoff);
	//} catch(const char *e) {
	//	std::cerr << e << std::endl;
	//	return 1;
	//}

	return 0;
}



