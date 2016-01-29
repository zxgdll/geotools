/*
 * planeadjust.cpp
 *
 *  Created on: Jan 16, 2016
 *      Author: rob
 */

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>

#include <geos/triangulate/IncrementalDelaunayTriangulator.h>
#include <geos/triangulate/quadedge/QuadEdgeSubdivision.h>
#include <geos/triangulate/quadedge/Vertex.h>
#include <geos/triangulate/quadedge/QuadEdge.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/Envelope.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include "Util.hpp"
#include "Raster.hpp"

#define PI 3.14159265358979323846

typedef geos::geom::Envelope Envelope;
typedef geos::triangulate::quadedge::QuadEdgeSubdivision QuadEdgeSubdivision;
typedef geos::triangulate::IncrementalDelaunayTriangulator IncrementalDelaunayTriangulator;
typedef geos::triangulate::quadedge::Vertex Vertex;
typedef geos::triangulate::quadedge::QuadEdge QuadEdge;
typedef geos::geom::Coordinate Coordinate;
typedef geos::geom::GeometryCollection GeometryCollection;
typedef geos::geom::GeometryFactory GeometryFactory;
typedef geos::geom::Geometry Geometry;

float _min(float a, float b) {
	return a > b ? b : a;
}

float _max(float a, float b) {
	return a < b ? b : a;
}

float _sq(float a) {
	return a*a;
}

class Point {
public:
	float p_x, p_y;
	float p_diff;
	Envelope *p_envelope;
	Point(float x, float y) :
		p_x(x),
		p_y(y),
		p_diff(nan("")),
		p_envelope(NULL) {

	}
	Point(Point &p) :
		Point(p.p_x, p.p_y) {
	}
	const Envelope *envelope() {
		if(!p_envelope)
			p_envelope = new Envelope(p_x, p_y, p_x, p_y);
		return p_envelope;
	}
	bool equals(const Point &p) const {
		return p_x == p.p_x && p_y == p.p_y;
	}
	~Point() {
		delete p_envelope;
	}
};

/**
 * Generates a set of sample points, one for each valid pixel in the mask.
 */
void generateSamplePoints(Raster<int> &mask, std::list<Point*> &samples) {
	for(int r = 0; r < mask.rows(); ++r) {
		for(int c = 0; c < mask.cols(); ++c) {
			//std::cerr << c << "," << r << ": " << mask.get(c, r) << std::endl;
			if(mask.get(c, r) != 0)
				samples.push_back(new Point(mask.toX(c), mask.toY(r)));
		}
	}
}

class samplesort {
public:
	samplesort(Raster<float> &r) : r(r) {
	}
	bool operator()(const Point *a, const Point *b) const {
		int ar = r.toBlockRow(a->p_y);
		int br = r.toBlockRow(b->p_y);
		int ac = r.toBlockCol(a->p_x);
		int bc = r.toBlockCol(b->p_x);
		if(ar != br) {
			return ar < br;
		} else {
			return ac < bc;
		}
	}

private:
	Raster<float> &r;
};

/**
 * Computes the raster differences for each sample point.
 */
void computeSampleDifferences(Raster<float> &base, Raster<float> &adj, std::list<Point*> &samples) {
	const samplesort ssort(base);
	samples.sort(ssort);
	for(std::list<Point*>::iterator pt = samples.begin(); pt != samples.end(); ++pt) {
		Point *p = *pt;
		p->p_diff = base.get(p->p_x, p->p_y) - adj.get(p->p_x, p->p_y);
		std::cerr << p->p_diff << " ";
	}
	samples.push_back(new Point(adj.toX(0), adj.toY(0)));
	samples.push_back(new Point(adj.toX(adj.cols() - 1), adj.toY(0)));
	samples.push_back(new Point(adj.toX(adj.cols() - 1), adj.toY(adj.rows() - 1)));
	samples.push_back(new Point(adj.toX(0), adj.toY(adj.rows() - 1)));
}

/**
 */
void adjust(std::string &basefile, std::string &adjfile, std::string &maskfile, std::string &outfile, int numSamples) {

	if(basefile == outfile || adjfile == outfile || maskfile == outfile)
		throw "The output file must not be the same as any input file.";

	if(numSamples < 1)
		throw "Too few samples.";

	GDALAllRegister();

	Envelope env;
	std::list<Point*> samples;

	{
		GDALDataset *maskDS = (GDALDataset *) GDALOpen(maskfile.c_str(), GA_ReadOnly);
		if(maskDS == NULL)
			throw "Failed to open adjustment file.";
		Raster<int> mask(maskDS, 1);
		generateSamplePoints(mask, samples);
		env.init(mask.toX(0) - 1.0, mask.toY(0) - 1.0, mask.toX(mask.cols()) + 1.0, mask.toY(mask.rows()) + 1.0);
	}

	GDALDataset *baseDS = (GDALDataset *) GDALOpen(basefile.c_str(), GA_ReadOnly);
	if(baseDS == NULL)
		throw "Failed to open base file.";
	GDALDataset *adjDS = (GDALDataset *) GDALOpen(adjfile.c_str(), GA_ReadOnly);
	if(adjDS == NULL)
		throw "Failed to open adjustment file.";
	{
		GDALDataset *outDS = adjDS->GetDriver()->CreateCopy(outfile.c_str(), adjDS, 1, NULL, NULL, NULL);
		if(outDS == NULL)
			throw "Failed to copy adjustment file.";
	}
	adjDS = (GDALDataset *) GDALOpen(outfile.c_str(), GA_Update);
	if(adjDS == NULL)
		throw "Failed to open output file.";

	Raster<float> base(baseDS, 1, false);
	Raster<float> adj(adjDS, 1, true);

	// TODO: Choose random sample from samples.
	computeSampleDifferences(base, adj, samples);

	std::map<Vertex, float> amap;

	GeometryFactory gf;
	QuadEdgeSubdivision qe(env, 1.0);
	IncrementalDelaunayTriangulator dt(&qe);

	for(std::list<Point*>::iterator pt = samples.begin(); pt != samples.end(); ++pt) {
		std::cerr << (*pt)->p_x << ", " << (*pt)->p_y << std::endl;
		const Vertex v((*pt)->p_x, (*pt)->p_y);
		amap[v] = (*pt)->p_diff;
		dt.insertSite(v);
	}


	for(int r = 0; r < adj.rows(); ++r) {
		for(int c = 0; c < adj.cols(); ++c) {
			
			float x = adj.toX(c);
			float y = adj.toY(r);
			const Vertex vc(x, y);

			// Add vertex

			dt.insertSite(vc);

			// Get vertex polygon and neighbours; compute areas
			std::map<QuadEdge*, double> areas;

			QuadEdge *e = qe.locate(vc);
			// Reverse the edge if the vc is at the wrong end.
			if(e->dest().equals(vc))
				e = &(e->sym());
			QuadEdge *e0 = e;

			qe.getVoronoiDiagram(gf);

			std::auto_ptr<Geometry> newCell = qe.getVoronoiCellPolygon(e0, gf);

			do {
				std::auto_ptr<Geometry> g = qe.getVoronoiCellPolygon(&(e0->sym()), gf);
				std::cerr << "Poly: " << g.get()->getArea() << std::endl;
				areas[e0] = g.get()->getArea();
				e0 = &(e0->oNext());
			} while(e0 != e);

			float z = 0.0;

			for(std::map<QuadEdge*, double>::iterator it = areas.begin(); it != areas.end(); ++it) {
				std::auto_ptr<Geometry> g = qe.getVoronoiCellPolygon(&(it->first->sym()), gf);
				std::cerr << "Poly: " << g.get()->getArea() << std::endl;
				z += amap[it->first->orig()] * g.get()->getArea() / areas[it->first];
			} 

			adj.set(x, y, z);
		}
	}


}

void usage() {
	std::cout << "Usage: rastfit <base file> <file to adjust> <mask file> <output file> <samples>" << std::endl;
	std::cout << "	This program will adjust one raster's elevations to match another's " << std::endl;
	std::cout << "	using an interpolated (IDW) adjustment to all pixels. The mask is " << std::endl;
	std::cout << "  required to limit the placement of samples where elevations are " << std::endl;
	std::cout << "	known to be good." << std::endl;
}

int main(int argc, char **argv) {

 	try {

 		if(argc < 6)
 			throw "Too few arguments.";

 		std::string basefile = argv[1];
 		std::string adjfile = argv[2];
 		std::string maskfile = argv[3];
  		std::string outfile = argv[4];
 		int samples = atoi(argv[5]);

 		adjust(basefile, adjfile, maskfile, outfile, samples);

 	} catch(const char *e) {
 		std::cerr << e << std::endl;
 		usage();
 		return 1;
 	}

 	return 0;
 }



