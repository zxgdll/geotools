/*
 * treecrowns.cpp
 *
 *  Created on: May 3, 2016
 *      Author: rob
 */

#include <queue>
#include <map>
#include <iostream>

#include "Raster.hpp"
#include "Vector.hpp"

class Node {
public:
	std::map<int, Node*> children;
	std::map<int, Node*> parents;
	int level;
	int id; 
	int row, col;

	Node(int id, int level, int col, int row) {
		this->id = id;
		this->level = level;
		this->col = col;
		this->row = row;
	}

	void addChild(Node *child) {
		if(children.find(child->id) == children.end())
			children[child->id] = child;
	}

	void addParent(Node *parent) {
		if(parents.find(parent->id) == parents.end())
			parents[parent->id] = parent;
	}

	int depth() {
		int l = 0;
		for(auto it = parents.begin(); it != parents.end(); ++it)
			l = g_max(l, it->second->depth());
		return l + 1;
	}

	~Node() {

	}
};


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
			weights[r * size + c] = (1 / (2 * G_PI * sigma * sigma)) * pow(E, -((x*x + y*y) / (2.0*sigma*sigma)));
		}
	}
}

/**
 * Smooth the raster and write the smoothed version to the output raster.
 */
// TODO: No accounting for nodata.
// TODO: Move to Raster.
void smooth(Grid<float> &raster, Grid<float> &smoothed, double sigma, int size) {
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

/** "Segment" the raster. By binning the heights into levels (starting with
 * 1). The range in heights is from minHeight to the raster's max. Anything 
 * below minHeight is set to 0.
 */
void segment(Grid<float> &src, Grid<int> &levels, Grid<int> &ids, int numLevels, float minHeight) {
	float min = minHeight;
	float max = src.max();
	float step = (max - min) / numLevels;
	levels.fill(0);
	for(unsigned int i = 0; i < src.size(); ++i) {
		if(src.get(i) < minHeight) {
			levels.set(i, 0);
		} else {
			levels.set(i, (int) ((src.get(i) - min) / step) + 1);
		}
	}

	// Set IDs to the negative value of the level patches.
	for(unsigned int i = 0; i < levels.size(); ++i)
		ids.set(i, -((int) levels.get(i)));

	// Flood fill the patches with new IDs -- a unique one for each patch.
	int id = 0;
	for(int r = 0; r < ids.rows(); ++r) {
		for(int c = 0; c < ids.cols(); ++c) {
			int value = ids.get(c, r);
			if(value < 0)
				ids.floodFill(c, r, value, ++id);
		}
	}
}

void buildGraph(Grid<int> &levels, Grid<int> &ids, std::map<int, Node*> &nodeMap) {
	for(int r = 0; r < levels.rows(); ++r) {
		for(int c = 0; c < levels.cols(); ++c) {
			
			int level = levels.get(c, r);
			int id = ids.get(c, r);
			Node *parent = nullptr;

			if(nodeMap.find(id) == nodeMap.end()) {
				parent = nodeMap[id] = new Node(id, level, c, r);
			} else {
				parent = nodeMap[id];
			}

			for(int rr = r - 1; rr < r + 2; ++rr) {
				for(int cc = c - 1; cc < c + 2; ++cc) {
					if((cc == c && rr == r) || cc < 0 || rr < 0 || 
						cc >= levels.cols() || rr >= levels.rows())
						continue;

					int id0 = ids.get(cc, rr);
					if(id0 == id)
						continue;

					int level0 = levels.get(cc, rr); // It's possible to have an equal level in the D8 case.
					Node *node = nullptr;

					if(nodeMap.find(id0) == nodeMap.end()) {
						node = nodeMap[id0] = new Node(id0, level0, cc, rr);
					} else {
						node = nodeMap[id0];
					}

					if(level0 < level) {
						parent->addChild(node);
						node->addParent(parent);
					} else if(level0 > level) {
						node->addChild(parent);
						parent->addParent(node);
					}
				}
			}
		}
	}
}

void printNode(Node *node, int level = 0, int depth = 1) {
	if(--depth < 0) return;
	char *pad = (char *) calloc(sizeof(char), level * 2 + 1);
	memset(pad, ' ', level * 2);
	std::cout << pad << "Node: " << node->id << ", " << node->level << "; " << node->col << ", " << node->row << std::endl;
	std::cout << pad << " - parents: " << node->parents.size() << std::endl;
	std::cout << pad << " - children: " << node->children.size() << std::endl;
	std::cout << pad << " -- parents: " << std::endl;
	for(auto it = node->parents.begin(); it != node->parents.end(); ++it) {
		printNode(it->second, level + 1, depth);
	}
	std::cout << pad << " -- chilren: " << std::endl;
	for(auto it = node->children.begin(); it != node->children.end(); ++it) {
		printNode(it->second, level + 1, depth);
	}
	std::cout << pad << "--------------------------------------" << std::endl;
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
	//if(crownshp.empty())
	//	throw "The crown output name cannot be empty.";
	//if(topshp.empty())
	//	throw "The treetop output name cannot be empty.";
	//if(smraster.empty())
	//	throw "The smoothed raster output name cannot be empty.";
	//if(smraster == topshp || topshp == crownshp || smraster == crownshp)
	//	throw "The output file names are in conflict.";

	std::cerr << "Loading " << inraster << std::endl;


	// Input raster.
	Raster<float> raster(inraster);
	MemRaster<float> smoothed(raster.cols(), raster.rows());
	MemRaster<int> levels(raster.cols(), raster.rows());
	levels.fill(0);
	levels.nodata(0);
	MemRaster<int> ids(raster.cols(), raster.rows());
	levels.fill(0);
	ids.nodata(0);
	
	// Smooth the input raster.
	std::cerr << "Smoothing..." << std::endl;
	smooth(raster, smoothed, 0.5, 3);

	// Segment the smoothed raster into patches.
	std::cerr << "Segmenting..." << std::endl;
	segment(smoothed, levels, ids, 10, 3.0);

	//Raster<float> smoothedr(std::string("/tmp/smoothed.tif"), raster);
	Raster<int> levelsr(std::string("/tmp/levels.tif"), raster);
	//Raster<int> idsr(std::string("/tmp/ids.tif"), raster);
	//smoothedr.writeBlock(smoothed);
	levelsr.writeBlock(levels);
	//idsr.writeBlock(ids);

	Raster<int> idsr1(std::string("/tmp/ids1.tif"), raster);
	idsr1.writeBlock(ids);

	std::cerr << "Building graph..." << std::endl;
	// Build a graph on the patches.
	std::map<int, Node*> nodes;
	buildGraph(levels, ids, nodes);

	ids.fill(0);

	std::vector<Node*> rootNodes;
	for(auto it = nodes.begin(); it != nodes.end(); ++it) {
		if(it->second->depth() == 2)
			rootNodes.push_back(it->second);
	}
	for(Node *node:rootNodes) {
		//printNode(node, 0, 2);
		levels.floodFill(node->col, node->row, node->level, 0, &ids, node->id);
	}

	Raster<int> idsr(std::string("/tmp/ids2.tif"), raster);
	idsr.writeBlock(ids);

	/*
	// Treetops vector file.
	remove(topshp.c_str());
	// Configure attributes.
	std::map<std::string, int> attribs;
	attribs["id"] = Vector::INTEGER;
	attribs["height"] = Vector::DOUBLE;
	Vector tops(topshp, Vector::POINT, proj, attribs);

	// Tree crown objects.
	std::vector<Crown> crowns;


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
	*/
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



