/*
 * rastfit
 *
 * This program adjusts a raster's elevations to match another's using natural
 * neighbour interpolation. It does this by collecting a configurable number of
 * random samples (constrained by an optional mask) and using these as sites in the 
 * triangulation. 
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

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/ch_jarvis.h>

#include "Util.hpp"
#include "Raster.hpp"
//#include "ShapeWriter.hpp"

typedef CGAL::Simple_cartesian<float>                                                 K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K>                  Vb; // Vertex can store its area
typedef CGAL::Triangulation_data_structure_2<Vb>                                      Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                                        Delaunay;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<Delaunay>                    AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<Delaunay>    AP;
typedef CGAL::Voronoi_diagram_2<Delaunay,AT,AP>                                       Voronoi;
typedef CGAL::Polygon_2<K>                                                            Polygon_2;
typedef CGAL::Iso_rectangle_2<K>                                                      Iso_rectangle_2;
typedef CGAL::Direction_2<K>                                                          Direction_2;

typedef Delaunay::Vertex_handle             DVertex_handle;
typedef Delaunay::Edge                      DEdge;
typedef Voronoi::Face_handle                VFace_handle;
typedef Voronoi::Face                       VFace;
typedef Voronoi::Ccb_halfedge_circulator    VCcb_halfedge_circulator;
typedef Voronoi::Locate_result              VLocate_result;
typedef Voronoi::Halfedge                   VHalfedge;
typedef Voronoi::Halfedge_handle            VHalfedge_handle;
typedef Voronoi::Vertex_handle              VVertex_handle;
typedef K::Point_2                          Point_2;
typedef K::Ray_2                            Ray_2;
typedef K::Segment_2                        Segment_2;

// For writing shapes for debugging.
//static ShapeWriter _sw;

float _abs(float v) {
	return v < 0.0 ? -v : v;
}

class Point {
public:
	float x, y, diff;
	Point() : x(nan("")), y(nan("")), diff(nan("")) {}
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
 * Sorts points so that they're in row/col order, 
 * which improves the read efficiency on the raster.
 */
class samplesort {
public:
	samplesort(Raster<float> &r) : r(r) {
	}
	bool operator()(std::unique_ptr<Point> &a, std::unique_ptr<Point> &b) const {
		int ar = r.toBlockRow(a->y);
		int br = r.toBlockRow(b->y);
		int ac = r.toBlockCol(a->x);
		int bc = r.toBlockCol(b->x);
		if(ar == br) {
			return ac < bc;
		} else {
			return ar < br;
		}
	}
private:
	Raster<float> &r;
};

/**
 * Shuffle and truncate the list of points, adding them to samples.
 */
void shuffle(std::vector<std::unique_ptr<Point> > &samp, std::list<std::unique_ptr<Point> > &samples, int numSamples) {
	std::random_shuffle(samp.begin(), samp.end());
	for(int i = 0; i < numSamples; ++i)
		samples.push_back(std::move(samp[i])); // TODO: Why?
	//samples.assign(samp.begin(), samp.begin() + numSamples);
}

/**
 * Computes the raster differences for each sample point.
 */
void computeSampleDifferences(Raster<float> &base, Raster<float> &adj, std::list<std::unique_ptr<Point> > &samples) {
	const samplesort ssort(base);
	samples.sort(ssort); // Sort because we want to read the raster efficiently
	for(auto pt = samples.begin(); pt != samples.end(); ++pt) {
		float a = base.get((*pt)->x, (*pt)->y);
		float b = adj.get((*pt)->x, (*pt)->y);
		(*pt)->diff = a - b;
	}
}

/**
 * Generates a set of sample points, limited to valid pixels in the mask.
 */
void generateMaskSamples(std::list<std::unique_ptr<Point> > &samples, Raster<char> &mask, int numSamples) {
	std::vector<std::unique_ptr<Point> > pts;
	for(int r = 0; r < mask.rows(); ++r) {
		for(int c = 0; c < mask.cols(); ++c) {
			if((int) mask.get(c, r) != 0) 
				pts.push_back(std::unique_ptr<Point>(new Point(mask.toX(c), mask.toY(r))));
		}
	}
	shuffle(pts, samples, numSamples);
}

/**
 * Generates a list of random samples within the bounds of a raster.
 */
void generateRandomSamples(std::list<std::unique_ptr<Point> > &samples, Raster<float> &rast, int numSamples) {
	do {
		float x = rast.toX((int) random() * rast.cols());
		float y = rast.toY((int) random() * rast.rows());
		samples.push_back(std::unique_ptr<Point>(new Point(x, y)));
	} while(--numSamples);
}

/**
 * Get a segment that represents the halfedge. If the halfedge is
 * finite, return an equivalent segment. If it is not, return a segment
 * cropped to the bounding rectangle.
 */
std::unique_ptr<Segment_2> getSegment(const VHalfedge &he, const Iso_rectangle_2 &boundary) {

	// TODO: Not using is_ray because of a bug. May not need this anymore.
	if(he.has_target() && he.has_source()) {

		Segment_2 s(he.source()->point(), he.target()->point());
		CGAL::Object inter;
		CGAL::Point_2<K> pi;
		CGAL::Segment_2<K> si;
		// The segment still might intersect the bounding rectangle. If it does,
		// return the clipped segment.
		for(int i = 0; i < 5; ++i) {
			// Rectangle side.
			Segment_2 s2(boundary.vertex(i), boundary.vertex((i + 1) % 4));
			if((inter = CGAL::intersection(s, s2))) {
				if(CGAL::assign(pi, inter)) {
					// If the intersection is a point, make a segment starting with the
					// point that is inside the boundary. If they're both outside, return null.
					if(boundary.bounded_side(s.source()) == CGAL::ON_BOUNDED_SIDE) {
						std::unique_ptr<Segment_2> seg(new Segment_2(s.source(), pi));
						return seg;
					} else if(boundary.bounded_side(s.target()) == CGAL::ON_BOUNDED_SIDE) {
						std::unique_ptr<Segment_2> seg(new Segment_2(pi, s.target()));
						return seg;
					} else {
						return NULL;
					}
				} else if(CGAL::assign(si, inter)) {
					// If the intersection is a segment, just return it.
					std::unique_ptr<Segment_2> seg(new Segment_2(si));
					return seg;			
				}
			}
		}
		// If the segment is inside the rectangle, return it.
		if(boundary.bounded_side(s.source()) == CGAL::ON_BOUNDED_SIDE) {
			std::unique_ptr<Segment_2> seg(new Segment_2(s));
			return seg;			
		}

	} else if(he.has_target() || he.has_source()) {

		// Get the dual of the halfedge and find its midpoint. This
		// helps determine the orientation of the ray.
		DEdge e = he.dual();
		DVertex_handle v1 = e.first->vertex((e.second + 1) % 3);
		DVertex_handle v2 = e.first->vertex((e.second + 2) % 3);

		// Direction perpendicular and to the right of v1, v2.
		Direction_2 dir(v1->point().y() - v2->point().y(), v2->point().x() - v1->point().x());
		// If the ray has a target, reverse the dir.
		Ray_2 ray(he.has_source() ? he.source()->point() : he.target()->point(), he.has_source() ? dir : -dir);

		// Check and return the intersection.
		CGAL::Object inter = CGAL::intersection(boundary, ray);
		CGAL::Point_2<K> p;
		CGAL::Segment_2<K> s;
		if(CGAL::assign(p, inter)) {
			std::unique_ptr<Segment_2> seg(new Segment_2(ray.source(), p));
			return seg;
		} else if(CGAL::assign(s, inter)) {
			std::unique_ptr<Segment_2> seg(new Segment_2(s));
			return seg;			
		}
	}

	return NULL;
}

/**
 * Returns the area of the clipped face.
 */
float faceArea(VFace f, Iso_rectangle_2 &boundary, Voronoi &vor) {
	std::list<Point_2> pts;
	VCcb_halfedge_circulator hc = f.ccb(), done(hc);
	do {
		std::unique_ptr<Segment_2> seg = getSegment(*hc, boundary);
		if(seg) {
			pts.push_back(seg.get()->source());
			pts.push_back(seg.get()->target());
		}
	} while(++hc != done);

	// Add the boundary vertices that are inside the current face.
	for(int i = 0; i < 4; ++i) {
		VLocate_result lr = vor.locate(boundary.vertex(i));
		VFace_handle *fh = boost::get<VFace_handle>(&lr);
		if(fh && **fh == f)
			pts.push_back(boundary.vertex(i));
	}

	// Build a convex hull of the points, and use the polygon to check the area.
	std::vector<Point_2> hull;
	CGAL::ch_jarvis(pts.begin(), pts.end(), std::back_inserter(hull));
	Polygon_2 poly(hull.begin(), hull.end());
	
	//_sw.put(poly);

	return _abs(poly.area());

}

/**
 * Compute an adjustment for adjfile, based on basefile, maskfile and a number of samples. Write the 
 * result to outfile.
 */
void adjust(std::string &basefile, std::string &adjfile, std::string &maskfile, std::string &outfile, int numSamples) {

	if(basefile == outfile || adjfile == outfile || maskfile == outfile)
		throw "The output file must not be the same as any input file.";

	if(numSamples < 1)
		throw "Too few samples.";

	// Initializes source, destination rasters.
	Raster<float> base(basefile, 1, false);
	Raster<float> adj(adjfile, 1, false);
	Raster<float> out = adj.copy(outfile, 1, true);

	std::list<std::unique_ptr<Point> > samples;

	// Generate a list of samples, shuffle it and clip it to the
	// desired count.
	if(maskfile.empty()) {
		generateRandomSamples(samples, adj, numSamples);
	} else {
		Raster<char> mask(maskfile, 1);
		generateMaskSamples(samples, mask, numSamples);
	}

	// Compute sample differences.
	computeSampleDifferences(base, adj, samples);

	// Build a boundary for clipping the voronoi.
	Iso_rectangle_2 boundary(
		Point_2(adj.toX(0) - 1000.0, adj.toY(0) + 1000.0),
		Point_2(adj.toX(adj.cols()) + 1000.0, adj.toY(adj.rows()) - 1000.0)
	);

	// Debug boundary.
	//_sw.put(boundary);

	// Start a delaunay triangulation.
	Delaunay dt;

	// Maps for site differences and face areas (by vertex ID).
	std::map<unsigned int, float> diffs;
	std::map<unsigned int, float> areas;

	// Build the delaunay triangulation on the samples, and associate the differences
	// with the vertices.
	unsigned int id = 0;
	for(auto pt = samples.begin(); pt != samples.end(); ++pt) {
		Point_2 p((*pt)->x, (*pt)->y);
		DVertex_handle h = dt.insert(p);
		h->info() = id;
		diffs[id] = (*pt)->diff;
		++id;
		// Debug point.
		//_sw.put(p);
	}

	// Pre-compute the areas of the original faces.
	Voronoi vt(dt);
	Voronoi::Face_iterator bf = vt.faces_begin();
	do {
		DVertex_handle h = bf->dual();
		float area = faceArea(*bf, boundary, vt);
		areas[h->info()] = area;
	} while(++bf != vt.faces_end());


	// Add the point for the initial cell to the triangulation.
	Point_2 vc(adj.toX(0), adj.toY(0));
	DVertex_handle vh = dt.insert(vc);

	for(int r = 0; r < adj.rows(); ++r) {
		std::cerr << "Row " << r << " of " << adj.rows() << std::endl;
		for(int c = 0; c < adj.cols(); ++c) {

			Point_2 vc(adj.toX(c), adj.toY(r));
			dt.move_if_no_collision(vh, vc);

			Voronoi vt0(dt);
			VLocate_result lr = vt0.locate(vh->point());
			VFace_handle fh = boost::get<VFace_handle>(lr);

			float area = faceArea(*fh, boundary, vt0);
			float z = adj.get(c, r), z0 = z;

			VCcb_halfedge_circulator hc = vt0.ccb_halfedges(fh), done(hc);
			do {
				VFace_handle f = hc->twin()->face();
				DVertex_handle v = f->dual();
				float a = faceArea(*f, boundary, vt0);
				unsigned int id = v->info();
				z0 += (areas[id] - a) / area * diffs[id];
			} while(++hc != done);

			out.set(c, r, z0);
		}
	}

	// Flush the debug shapes.
	//_sw.write();

}

void usage() {
	std::cout << "Usage: rastfit <base file> <file to adjust> <mask file | \"-\"> <output file> <samples>" << std::endl;
	std::cout << "	This program will adjust one raster's elevations to match another's " << std::endl;
	std::cout << "	using a natural neighbours interpoled adjustment to all pixels. The mask is " << std::endl;
	std::cout << "  optional and used to limit the placement of samples to where elevations are " << std::endl;
	std::cout << "	known to be good." << std::endl;
}

int main(int argc, char **argv) {

 	try {

 		if(argc < 6)
 			throw "Too few arguments.";

 		//_sw.off();

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



