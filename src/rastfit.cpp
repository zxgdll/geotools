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
#include <algorithm>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include "Util.hpp"
#include "Raster.hpp"
#include "ShapeWriter.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel 							K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K>					Vb;			// Vertex can store its area
typedef CGAL::Triangulation_data_structure_2<Vb>										Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> 											Delaunay;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<Delaunay> 						AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<Delaunay> 		AP;
typedef CGAL::Voronoi_diagram_2<Delaunay,AT,AP> 										Voronoi;
typedef CGAL::Polygon_2<K>																Polygon_2;
typedef CGAL::Iso_rectangle_2<K>														Iso_rectangle_2;

typedef Delaunay::Edge_circulator 			DEdge_circulator;
typedef Delaunay::Vertex_handle 			DVertex_handle;
typedef Delaunay::Edge						DEdge;
typedef Voronoi::Face_handle 				VFace_handle;
typedef Voronoi::Face 						VFace;
typedef Voronoi::Ccb_halfedge_circulator	VCcb_halfedge_circulator;
typedef Voronoi::Locate_result				VLocate_result;
typedef Voronoi::Halfedge					VHalfedge;
typedef Voronoi::Halfedge_handle			VHalfedge_handle;
typedef Voronoi::Vertex_handle				VVertex_handle;
typedef K::Point_2 							Point_2;
typedef K::Ray_2							Ray_2;
typedef K::Segment_2						Segment_2;
typedef K::Line_2							Line_2;

#define PI 3.14159265358979323846

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
	float x, y, diff;
	Point(float x, float y) :
		x(x),
		y(y),
		diff(nan("")) {
	}
	Point(float x, float y, float z) :
		x(x),
		y(y),
		diff(z) {
	}
	Point(Point &p) :
		Point(p.x, p.y, p.diff) {
	}
	bool equals(const Point &p) const {
		return x == p.x && y == p.y;
	}
	~Point() {
	}
};

/**
 * Generates a set of sample points, one for each valid pixel in the mask.
 */
void generateSamplePoints(Raster<char> &mask, std::vector<Point*> &samples) {
	for(int r = 0; r < mask.rows(); ++r) {
		for(int c = 0; c < mask.cols(); ++c) {
			if((int) mask.get(c, r) != 0)
				samples.push_back(new Point(mask.toX(c), mask.toY(r)));
		}
	}
}

class samplesort {
public:
	samplesort(Raster<float> &r) : r(r) {
	}
	bool operator()(const Point *a, const Point *b) const {
		int ar = r.toBlockRow(a->y);
		int br = r.toBlockRow(b->y);
		int ac = r.toBlockCol(a->x);
		int bc = r.toBlockCol(b->x);
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
	samples.sort(ssort); // Sort because we want to read the raster efficiently
	for(std::list<Point*>::iterator pt = samples.begin(); pt != samples.end(); ++pt) {
		Point *p = *pt;
		float a = base.get(p->x, p->y);
		float b = adj.get(p->x, p->y);
		p->diff = a - b;
	}
}

float _abs(float v) {
	return v < 0.0 ? -v : v;
}

static ShapeWriter _sw;

/**
 * Get a segment that represents the halfedge. If the halfedge is
 * finite, return an equivalent segment. If it is not, return a segment
 * cropped to the bounding rectangle.
 */
std::unique_ptr<Segment_2> getSegment(const VHalfedge &he, const Iso_rectangle_2 &boundary, const Voronoi &vor) {
	CGAL::Object obj;
	if(!he.is_unbounded()) {
		Segment_2 seg(he.source()->point(), he.target()->point());
		obj = CGAL::intersection(boundary, seg);
	} else if(he.has_source()) {
		// Get the dual of the halfedge and find its midpoint. This
		// helps determine the orientation of the ray.
		DEdge e = he.dual();
		DVertex_handle v1 = e.first->vertex((e.second + 2) % 3);
		DVertex_handle v2 = e.first->vertex((e.second + 1) % 3);
		float x1 = v1->point().x(), x2 = v2->point().x();
		float y1 = v1->point().y(), y2 = v2->point().y();
		Point_2 mid((x1 + x2) / 2.0, (y1 + y2) / 2.0);

		// If the bisector does not cross the ray, 
		// if it does, reverse the ray. TODO: This is broken logic.
		bool reverse = false;
		VLocate_result lr = vor.locate(mid);
		VHalfedge_handle *hh = boost::get<VHalfedge_handle>(&lr);	
		if(!hh) {
			reverse = true;
		} else if(**hh != he) {
			// TODO: This should probably be active...
			//reverse = true;
		}
		
		Point_2 pt = he.source()->point();
		if(reverse) {
			// If we're reversing the ray, move the midpoint here.
			float x0 = mid.x() + (pt.x() - mid.x()) * 2;
			float y0 = mid.y() + (pt.y() - mid.y()) * 2;
			mid = Point_2(x0, y0);
		}
		Ray_2 r(pt, mid);
		obj = CGAL::intersection(boundary, r);
	}

	if(obj == NULL) // The ray doesn't intersect the boundary.
		return NULL;

	// Obj could be a segment or a point. If it's a point, 
	// use the source of the ray.
	if(const Segment_2 *se = CGAL::object_cast<Segment_2>(&obj)) {
		std::unique_ptr<Segment_2> ptp(new Segment_2(*se));
		//_sw.put(*ptp.get());
		return ptp;
	} else if(const Point_2 *p = CGAL::object_cast<Point_2>(&obj)) {
		std::unique_ptr<Segment_2> ptp(new Segment_2(he.source()->point(), *p));
		//_sw.put(*ptp.get());
		return ptp;
	}

	return NULL;

}

/**
 * Returns the area of a face.
 */
float faceArea(VFace f, Iso_rectangle_2 &boundary, Voronoi &vor) {
	std::list<Point_2> pts;
	VCcb_halfedge_circulator hc = f.ccb(), done(hc);
	do {
		std::unique_ptr<Segment_2> seg = getSegment(*hc, boundary, vor);
		if(seg) {
			pts.push_back(seg->source());
			pts.push_back(seg->target());
			_sw.put(*seg.get());
		}
	} while(++hc != done);
	std::unique_ptr<Segment_2> seg = getSegment(*done, boundary, vor);
	if(seg) {
		pts.push_back(seg->source());
		pts.push_back(seg->target());
		_sw.put(*seg.get());
	}
	Polygon_2 poly(pts.begin(), pts.end());
	_sw.put(poly);
	return _abs(poly.area());
}

/**
 */
void adjust(std::string &basefile, std::string &adjfile, std::string &maskfile, std::string &outfile, int numSamples) {

	if(basefile == outfile || adjfile == outfile || maskfile == outfile)
		throw "The output file must not be the same as any input file.";

	if(numSamples < 1)
		throw "Too few samples.";

	std::list<Point*> samples;

	{
		Raster<char> mask(maskfile, 1);
		std::vector<Point*> samp;
		generateSamplePoints(mask, samp);
		std::random_shuffle(samp.begin(), samp.end());
		int i = 0;
		for(std::vector<Point*>::iterator it = samp.begin(); it != samp.end(); ++it) {
			if(i < numSamples) {
				samples.push_back(*it);
			} else {
				delete *it;
			}
			++i;
		}
	}

	Raster<float> base(basefile, 1, false);
	Raster<float> adj(adjfile, 1, false);
	Raster<float> out(outfile, 1, true);

	// TODO: Choose random sample from samples.
	computeSampleDifferences(base, adj, samples);

	Iso_rectangle_2 boundary(
		Point_2(adj.toX(0) - 1000.0, adj.toY(0) + 1000.0),
		Point_2(adj.toX(adj.cols()) + 1000.0, adj.toY(adj.rows()) - 1000.0)
	);

	_sw.put(boundary);

	Delaunay dt;

	std::map<unsigned int, float> diffs;
	std::map<unsigned int, float> areas;

	unsigned int id = 0;
	for(std::list<Point*>::iterator pt = samples.begin(); pt != samples.end(); ++pt) {
		Point_2 p((*pt)->x, (*pt)->y);

		_sw.put(p);

		DVertex_handle h = dt.insert(p);
		h->info() = id;
		diffs[id] = (*pt)->diff;
		std::cerr << "diff " << id << ": " << diffs[id] << std::endl;
		++id;
	}

	Voronoi vt(dt);
	Voronoi::Face_iterator bf = vt.faces_begin();
	do {
		DVertex_handle h = bf->dual();
		//if(dt.is_infinite(h))
		//	continue;
		float area = faceArea(*bf, boundary, vt);
		areas[h->info()] = area;
		std::cerr << "area " << h->info() << ": " << area << std::endl;
	} while(++bf != vt.faces_end());

	Point_2 vc(adj.toX(0), adj.toY(0));
	DVertex_handle vh = dt.insert(vc);

	_sw.write();

	for(int r = 0; r < adj.rows(); ++r) {
		std::cerr << r << std::endl;
		for(int c = 0; c < adj.cols(); ++c) {

			Point_2 vc(adj.toX(c), adj.toY(r));
			dt.move(vh, vc);

			Voronoi vt0(dt);
			VLocate_result lr = vt0.locate(vh->point());
			VFace_handle fh = boost::get<VFace_handle>(lr);

			float area = faceArea(*fh, boundary, vt0);
			float z = adj.get(c, r);
			std::cerr << "z1: " << z << ", " << area << std::endl;

			VCcb_halfedge_circulator hc = vt0.ccb_halfedges(fh), done(hc);
			float x = 0.0;
			do {
				VFace_handle f = hc->twin()->face();
				DVertex_handle v = f->dual();
				float a = faceArea(*f, boundary, vt0);
				unsigned int id = v->info();
				if(a > 0.0) {
					std::cerr << "area " << id << ": " << areas[id] << ", " << a << ", " << (areas[id] - a) << ", " << area << ", " << ((areas[id] - a) / area) << ", " << diffs[id] << std::endl;
					x += (areas[id] - a) / area;
					z -= (areas[id] - a) / area * diffs[id];
				}
			} while(++hc != done);

			std::cerr << "z2: " << z << ", " << x << std::endl;
			out.set(c, r, z);
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



