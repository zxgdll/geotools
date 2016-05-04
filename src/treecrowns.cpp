/*
 * treecrowns.cpp
 *
 *  Created on: May 3, 2016
 *      Author: rob
 */

#include <queue>

#include "Raster.hpp"

// Create ID grid (int)
// Find maxima
// Iterate over maxima;
//   create ID, add to queue
// For each in queue:
//   check surrounding pixels
//     if higher, stop
//     if lower than [dropoff], stop
//     else set id, add to queue

static unsigned int _crown_id = 0;

class Px {
public:
	int col;
	int row;
	unsigned int id;
	Px(int col, int row, unsigned int id) {
		this->col = col;
		this->row = row;
		this->id = id;
	}
};

class Crown {
public:
	unsigned int id;
	int col;
	int row;
	float height;
	Crown(int col, int row, float height) {
		this->col = col;
		this->row = row;
		this->height = height;
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

bool isMaximum(Raster<float> &raster, int col, int row, float *height) {
	static int ks = 7;
	static int ks0 = (int) ks / 2;
	if(col < ks0 || row < ks0 || col >= raster.cols() - ks0 || row >= raster.rows() - ks0)
		return false;
	static Grid<float> k(ks, ks);
	raster.readBlock(col - ks0, row - ks0, ks, ks, k);
	*height = k.get(ks0, ks0);
	for(int r = 0; r < ks; ++r) {
		for(int c = 0; c < ks; ++c) {
			if(r == ks0 && c == ks0) continue;
			if(k.get(c, r) > *height)
				return false;
		}
	}
	return true;
}

void findMaxima(Raster<float> &raster, std::vector<Crown> &crowns) {
	float height = nan("");
	for(int r = 0; r < raster.rows(); ++r) {
		for(int c = 0; c < raster.cols(); ++c) {
			if(isMaximum(raster, c, r, &height))
				crowns.push_back(Crown(c, r, height));
		}
	}
}

void treetops(std::string &inraster, std::string &outshp, double dropoff) {

	static Grid<float> k(3, 3);
	Raster<float> raster(inraster);
	Grid<int> ids(raster.cols(), raster.rows());
	ids.fill(0);

	std::vector<Crown> crowns;
	findMaxima(raster, crowns);

	//std::sort(crowns.begin(), crowns.end(), crownReverseSort);

	std::queue<Px> cells;
	for(Crown &crn:crowns)
		cells.push(Px(crn.col, crn.row, crn.id));

	while(!cells.empty()) {
		Px px = cells.front();
		cells.pop();
		if(px.col < 1 || px.row < 1 || px.col >= raster.cols() - 1 || px.row >= raster.rows() - 1)
			continue;
		if(ids.get(px.col, px.row) != 0)
			continue;
		raster.readBlock(px.col - 1, px.row - 1, 3, 3, k);
		float height = k.get(1, 1);
		ids.set(px.col, px.row, px.id);
		for(int r = 0; r < 3; ++r) {
			for(int c = 0; c < 3; ++c) {
				if(c == 1 && r ==1) continue;
				float height0 = k.get(c, r);
				if(height0 > height) continue;
				if(height0 < (height * dropoff)) continue;
				//ids.set(px.col - 1 + c, px.row - 1 + r, px.id);
				cells.push(Px(px.col - 1 + c, px.row - 1 + r, px.id));
			}
		}
	}

	std::string proj;
	raster.projection(proj);
	Raster<int> output(outshp, raster.minx(), raster.miny(), raster.maxx(), raster.maxy(),
				raster.resolution(), 0, proj);
	output.writeBlock(ids);

}

int main(int argc, char **argv) {

	std::string inraster;
	std::string outshp;
	double dropoff = 0.65; // Percentage of elevation to drop without stop.

	int i = 1;
	for(; i < argc; ++i) {
		std::string arg = argv[i];
		if(arg == "-d") {
			dropoff = atof(argv[++i]);
		} else if(arg == "-i") {
			inraster = argv[++i];
		} else if(arg == "-o") {
			outshp = argv[++i];
		}
	}

	treetops(inraster, outshp, dropoff);

}



