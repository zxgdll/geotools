/*
 * rastfit
 *
 * This program adjusts a raster's elevations to match another's using natural
 * neighbour interpolation or plane fitting. It does this by collecting a configurable
 * number of random samples constrained by an optional mask and using these as sites
 * to compute differences.
 *
 *  Created on: Jan 16, 2016
 *  Author: rob
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

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/ch_jarvis.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/number_utils.h>

#include "Util.hpp"
#include "Raster.hpp"
//#include "ShapeWriter.hpp"

typedef CGAL::Exact_predicates_exact_constructions_kernel                             K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K>                  Vb; // Vertex can store its area
typedef CGAL::Triangulation_data_structure_2<Vb>                                      Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                                        Delaunay;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<Delaunay>                    AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<Delaunay>    AP;
typedef CGAL::Voronoi_diagram_2<Delaunay,AT,AP>                                       Voronoi;
typedef CGAL::Polygon_2<K>                                                            Polygon_2;
typedef CGAL::Polygon_with_holes_2<K>                                                 Polygon_with_holes_2;
typedef CGAL::Iso_rectangle_2<K>                                                      Iso_rectangle_2;
typedef CGAL::Direction_2<K>                                                          Direction_2;
typedef CGAL::Point_2<K>                          Point_2;
typedef CGAL::Ray_2<K>                            Ray_2;
typedef CGAL::Segment_2<K>                        Segment_2;

typedef Delaunay::Vertex_handle             DVertex_handle;
typedef Delaunay::Edge                      DEdge;
typedef Voronoi::Face_handle                VFace_handle;
typedef Voronoi::Face                       VFace;
typedef Voronoi::Ccb_halfedge_circulator    VCcb_halfedge_circulator;
typedef Voronoi::Locate_result              VLocate_result;
typedef Voronoi::Halfedge                   VHalfedge;
typedef Voronoi::Halfedge_handle            VHalfedge_handle;
typedef Voronoi::Vertex_handle              VVertex_handle;

//ShapeWriter sw("/home/rob/Documents/geotools/data/out.shp");

double _abs(double v) {
	return v < 0.0 ? -v : v;
}

double _min(double a, double b) {
	return a > b ? b : a;
}

double _max(double a, double b) {
	return a < b ? b : a;
}

double _random() {
	double r = ((double) std::rand()) / RAND_MAX;
	return r;
}

class Point {
public:
	double x, y, diff;
	Point() : x(nan("")), y(nan("")), diff(nan("")) {}
	Point(double x, double y) :
		x(x),
		y(y),
		diff(nan("")) {
	}
	Point(double x, double y, double z) :
		x(x),
		y(y),
		diff(z) {
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
		return (a.get(), b.get());
	}
	bool operator()(Point &a, Point &b) const {
		int ar = r.toBlockRow(a.y);
		int br = r.toBlockRow(b.y);
		int ac = r.toBlockCol(a.x);
		int bc = r.toBlockCol(b.x);
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
 * Computes the raster differences for each sample point.
 */
void computeSampleDifferences(Raster<float> &base, Raster<float> &adj, std::list<Point> &samples) {
	const samplesort ssort(base);
	samples.sort(ssort); // Sort because we want to read the raster efficiently
	for(auto pt = samples.begin(); pt != samples.end(); ++pt) {
		float a = base.get(pt->x, pt->y);
		float b = adj.get(pt->x, pt->y);
		if(a != base.nodata() && b != adj.nodata()) {
			pt->diff = a - b;
		} else {
			pt->diff = adj.nodata();
		}
	}
}

/**
 * Generates a set of sample points, limited to valid pixels in the mask.
 */
void generateMaskSamples(std::list<Point> &samples, Raster<char> &mask, 
	Raster<float> &base, Raster<float> &adj, unsigned int numSamples) {
	std::vector<Point> pts;
	for(int r = 0; r < mask.rows(); ++r) {
		for(int c = 0; c < mask.cols(); ++c) {
			if(mask.get(c, r) != 0 && adj.isValid(mask.toX(c), mask.toY(r))
					&& base.isValid(mask.toX(c), mask.toY(r))) {
				pts.push_back(Point(mask.toX(c), mask.toY(r)));
			}
		}
	}
	std::cerr << "Selecting from " << pts.size() << " points" << std::endl;
	std::random_shuffle(pts.begin(), pts.end());
	samples.assign(pts.begin(), pts.begin() + (numSamples > pts.size() ? pts.size() : numSamples));
}

/**
 * Generates a list of random samples within the bounds of a raster.
 */
void generateRandomSamples(std::list<Point> &samples, Raster<float> &a, Raster<float> &b, 
	double *bounds, unsigned int numSamples) {
	while(numSamples-- > 0) {
		double x = bounds[0] + (bounds[2] - bounds[0]) * _random();
		double y = bounds[1] + (bounds[3] - bounds[1]) * _random();
		if(a.isValid(x, y) && b.isValid(x, y)) {
			Point p(x, y);
			samples.push_back(p);
		}
	}
}

/**
 * Get or generate the output points of the halfedge. If the edge is unbounded,
 * make a bounded segment using the d argument, which should be a reasonable distance
 * outside the bounds of the region of interest.
 */
void getSegmentPoints(const VHalfedge &he, std::list<Point_2> &pts, double d) {
	if(!he.is_ray()) {
		// If it's not a ray, it can just be returned as a segment.
		pts.push_back(he.source()->point());
		pts.push_back(he.target()->point());
	} else {
		// Get the dual of the halfedge and find its midpoint. This
		// helps determine the orientation of the ray.
		DEdge e = he.dual();
		DVertex_handle v1 = e.first->vertex((e.second + 1) % 3);
		DVertex_handle v2 = e.first->vertex((e.second + 2) % 3);
		// Direction perpendicular and to the right of v1, v2.
		Direction_2 dir(v1->point().y() - v2->point().y(), v2->point().x() - v1->point().x());
		// If the ray has a target, reverse the dir.
		if(!he.has_source())
			dir = -dir;
		double a = atan2(CGAL::to_double(dir.dy()), CGAL::to_double(dir.dx()));
		Point_2 pt1 = he.has_source() ? he.source()->point() : he.target()->point();
		Point_2 pt2(pt1.x() + d * cos(a), pt1.y() + d * sin(a));
		pts.push_back(pt1);
		pts.push_back(pt2);
	}
}

/**
 * Convert a Voronoi face to a polygon clipped to the given boundary.
 */
double faceArea(const VFace &f, const Voronoi &vor, const Polygon_2 &bounds) {
	// Compute a length for unbounded segments.
	double d = _max(bounds.bbox().xmax() - bounds.bbox().xmin(),
			bounds.bbox().ymax() - bounds.bbox().ymin()) * 2.0;
	std::list<Point_2> pts;
	// Build segments from the face edges.
	VCcb_halfedge_circulator hc = f.ccb(), done(hc);
	do {
		getSegmentPoints(*hc, pts, d);
	} while(++hc != done);

	// Add the boundary vertices that are inside the current face.
	// TODO: This is a hack but it'll do for now.
	for(int i = 0; i < 4; ++i) {
		VLocate_result lr = vor.locate(bounds.vertex(i));
		VFace_handle *fh = boost::get<VFace_handle>(&lr);
		if(fh && **fh == f)
			pts.push_back(bounds.vertex(i));
	}

	// Build a convex hull of the points. All cropped voronoi cells are convex.
	std::vector<Point_2> hull;
	CGAL::ch_jarvis(pts.begin(), pts.end(), std::back_inserter(hull));
	Polygon_2 poly(hull.begin(), hull.end());

	// Get the polygon intersection. There can be zero or 1.
	std::vector<Polygon_with_holes_2> out;
	CGAL::intersection(poly, bounds, std::back_inserter(out));

	// Return the poly, or none.
	if(out.size() > 0) {
		return CGAL::to_double(out[0].outer_boundary().area());
	} else {
		return 0.0;
	}
}

void printSamples(std::list<Point> &samples) {
	std::cout << "x,y,diff" << std::endl;
	for(auto it = samples.begin(); it != samples.end(); ++it) {
		std::cout << std::setprecision(9) << it->x << "," << it->y << "," << it->diff << std::endl;
	}
}

void computeBounds(Raster<float> &a, Raster<float> &b, double *bounds) {
	bounds[0] = _max(a.toX(0), b.toX(0));
	bounds[1] = _max(a.toY(a.rows()), b.toY(b.rows()));
	bounds[2] = _min(a.toX(a.cols()), b.toX(b.cols()));
	bounds[3] = _min(a.toY(0), b.toY(0));
}

/**
 * Compute an adjustment for adjfile, based on basefile, maskfile and a number of samples. Write the 
 * result to outfile.
 */
void adjustNN(std::string &basefile, std::string &adjfile, std::string &maskfile,
		std::string &outfile, unsigned int numSamples, double resolution) {

	if(basefile == outfile || adjfile == outfile || maskfile == outfile)
		throw "The output file must not be the same as any input file.";

	if(numSamples < 2)
		throw "Too few samples.";

	if(resolution <= 0.0)
		throw "Invalid resolution.";

	// Initializes source, destination rasters.
	Raster<float> base(basefile, 1, false);
	Raster<float> adj(adjfile, 1, false);
	std::string proj;
	adj.projection(proj);
	Raster<float> out(outfile, adj.minx(), adj.miny(), adj.maxx(), adj.maxy(), resolution, proj.c_str());

	base.get(0,0);
	adj.get(0,0);
	std::list<Point> samples;

	// Generate a list of samples, shuffle it and clip it to the
	// desired count.
	if(maskfile.empty() || maskfile == "-") {
		double bounds[4];
		computeBounds(base, adj, bounds);
		generateRandomSamples(samples, adj, base, bounds, numSamples);
	} else {
		std::cerr << "Generating samples." << std::endl;
		Raster<char> mask(maskfile, 1);
		generateMaskSamples(samples, mask, adj, base, numSamples);
	}

	if(samples.size() == 0)
		throw "No samples collected.";

	std::cerr << "Computing sample differences." << std::endl;
	// Compute sample differences.
	computeSampleDifferences(base, adj, samples);

	printSamples(samples);

	// Build a boundary for clipping the voronoi.
	std::list<Point_2> pts;
	pts.push_back(Point_2(adj.maxx() + 1000.0, adj.maxy() + 1000.0));
	pts.push_back(Point_2(adj.maxx() + 1000.0, adj.miny() - 1000.0));
	pts.push_back(Point_2(adj.minx() - 1000.0, adj.miny() - 1000.0));
	pts.push_back(Point_2(adj.minx() - 1000.0, adj.maxy() + 1000.0));
	pts.reverse();
	Polygon_2 bounds(pts.begin(), pts.end());

	// Start a delaunay triangulation.
	Delaunay dt;

	// Maps for site differences and face areas (by vertex ID).
	std::map<unsigned int, double> diffs;
	std::map<unsigned int, double> areas;

	// Build the delaunay triangulation on the samples, and associate the differences
	// with the vertices (which are sample sites.)
	unsigned int id = 0;
	for(auto pt = samples.begin(); pt != samples.end(); ++pt) {
		Point_2 p(pt->x, pt->y);
		DVertex_handle h = dt.insert(p);
		h->info() = id;
		diffs[id] = pt->diff;
		++id;
	}


	// Pre-compute the areas of the original faces.
	Voronoi vt(dt);
	Voronoi::Face_iterator bf = vt.faces_begin();
	do {
		DVertex_handle h = bf->dual();
		double area = faceArea(*bf, vt, bounds);
		if(area <= 0.0)
			throw "Invalid area.";
		areas[h->info()] = area;
	} while(++bf != vt.faces_end());

	// Add the point for the initial cell to the triangulation.
	Point_2 vc(out.toX(0), out.toY(0));
	DVertex_handle vh = dt.insert(vc);

	for(int r = 0; r < out.rows(); ++r) {
		std::cerr << "Row " << r << " of " << out.rows() << std::endl;
		for(int c = 0; c < out.cols(); ++c) {

			Point_2 vc(out.toX(c), out.toY(r));
			dt.move_if_no_collision(vh, vc);

			Voronoi vt0(dt);

			VLocate_result lr = vt0.locate(vh->point());
			VFace_handle fh = boost::get<VFace_handle>(lr);

			double area = faceArea(*fh, vt0, bounds);
			if(area == 0.0)
				throw "Invalid face area.";
			double z0 = 0.0;
			VCcb_halfedge_circulator hc = vt0.ccb_halfedges(fh), done(hc);
			do {
				VFace_handle f = hc->twin()->face();
				DVertex_handle v = f->dual();
				double a = faceArea(*f, vt0, bounds);
				unsigned int id = v->info();
				if((areas[id] - a) < -0.00001) // TODO: Precision problem -- look into it.
					throw "The modified area cannot be larger than the original area.";
				z0 += ((areas[id] - a) / area) * diffs[id];
			} while(++hc != done);
			out.set(c, r, z0);
		}
	}
}

void determinant3(double *d, double **m) {
	*d = m[0][0] * m[1][1] * m[2][2]
		+ m[0][1] * m[1][2] * m[2][0]
		+ m[0][2] * m[1][0] * m[2][1]
		- m[0][2] * m[1][1] * m[2][0]
		- m[0][1] * m[1][0] * m[2][2]
		- m[0][0] * m[1][2] * m[2][1];
}

void scaleadjoint3(double **a, double s, double **m) {
   a[0][0] = (s) * (m[1][1] * m[2][2] - m[1][2] * m[2][1]);
   a[1][0] = (s) * (m[1][2] * m[2][0] - m[1][0] * m[2][2]);
   a[2][0] = (s) * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
   a[0][1] = (s) * (m[0][2] * m[2][1] - m[0][1] * m[2][2]);
   a[1][1] = (s) * (m[0][0] * m[2][2] - m[0][2] * m[2][0]);
   a[2][1] = (s) * (m[0][1] * m[2][0] - m[0][0] * m[2][1]);
   a[0][2] = (s) * (m[0][1] * m[1][2] - m[0][2] * m[1][1]);
   a[1][2] = (s) * (m[0][2] * m[1][0] - m[0][0] * m[1][2]);
   a[2][2] = (s) * (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
}
 
double minvert3(double **a, double **b) {
	double det = 0;
	determinant3(&det, a);
	std::cerr << "det " << det << std::endl;
	double tmp = 1.0 / det;
	scaleadjoint3(b, tmp, a);
	return det;
}

void mtranspose(double **src, int maj, int min, double **dst) {
	for(int r = 0; r < maj; ++r) {
		for(int c = 0; c < min; ++c) {
			dst[c][r] = src[r][c];
		}
	}
}

void mmult(double **a, int amaj, int amin, double **b, int bmaj, int bmin, double **dst) {
	if(amin != bmaj)
		throw "Incompatible row/column sizes.";
	for(int ar = 0; ar < amaj; ++ar) {
		for(int bc = 0; bc < bmin; ++bc) {
			dst[ar][bc] = 0.0;
			for(int i = 0; i < amin; ++i)
				dst[ar][bc] += a[ar][i] * b[i][bc];
		}
	}
}

void mcentroid(double **mtx, int rows, int xcol, int ycol, double *x, double *y) {
	*x = 0, *y = 0;
	for(int r = 0; r < rows; ++r) {
		*x += mtx[r][xcol];
		*y += mtx[r][ycol];
	}
	*x /= rows;
	*y /= rows;
	for(int r = 0; r < rows; ++r) {
		mtx[r][xcol] -= *x;
		mtx[r][ycol] -= *y;
	}
}

double **minit(int rows, int cols) {
	double **mtx = new double*[rows];
	for(int r = 0; r < rows; ++r)
		mtx[r] = new double[cols];
	return mtx;
}

void mfree(double **mtx, int rows) {
	for(int r = 0; r < rows; ++r)
		delete mtx[r];
	delete mtx;
}

void mprint(double **mtx, int rows, int cols) {
	std::cerr << "mtx " << rows << "," << cols << std::endl;
	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c) {
			std::cerr << mtx[r][c] << " ";
		}
		std::cerr << std::endl;
	}
}

void adjustPlanar(std::string &basefile, std::string &adjfile, std::string &maskfile,
		std::string &outfile, unsigned int numSamples, double resolution) {

	if(basefile == outfile || adjfile == outfile || maskfile == outfile)
		throw "The output file must not be the same as any input file.";

	if(numSamples < 2)
		throw "Too few samples.";

	if(resolution <= 0.0)
		throw "Invalid resolution.";

	// Initializes source, destination rasters.
	Raster<float> base(basefile, 1, false);
	Raster<float> adj(adjfile, 1, false);
	Raster<char> mask(maskfile, 1, false);
	std::string proj;
	adj.projection(proj);
	Raster<float> out(outfile, adj.minx(), adj.miny(), adj.maxx(), adj.maxy(), resolution, proj.c_str());

	base.get(0,0);
	adj.get(0,0);
	std::list<Point> samples;

	// Generate a list of samples, shuffle it and clip it to the
	// desired count.
	std::cerr << "Generating samples." << std::endl;
	if(maskfile.empty() || maskfile == "-") {
		double bounds[4];
		computeBounds(base, adj, bounds);
		generateRandomSamples(samples, adj, base, bounds, numSamples);
	} else {
		generateMaskSamples(samples, mask, adj, base, numSamples);
	}

	if(samples.size() == 0)
		throw "No samples collected.";

	std::cerr << "Computing sample differences." << std::endl;
	// Compute sample differences.
	computeSampleDifferences(base, adj, samples);

	//printSamples(samples);

	double **xymtx = minit(samples.size(), 3);
	double **zmtx = minit(samples.size(), 1);
	int i = 0;
	for(auto it = samples.begin(); it != samples.end(); ++it) {
		xymtx[i][0] = 1.0;
		xymtx[i][1] = it->x;
		xymtx[i][2] = it->y;
		zmtx[i][0] = it->diff;
		++i;
	}

	double cx, cy;
	mcentroid(xymtx, samples.size(), 1, 2, &cx, &cy);

	double **xymtxt = minit(3, samples.size());
	mtranspose(xymtx, samples.size(), 3, xymtxt);

	double **a = minit(3, 3);
	mmult(xymtxt, 3, samples.size(), xymtx, samples.size(), 3, a);

	double **ai = minit(3, 3);
	minvert3(a, ai);

	double **b = minit(3, 3);
	mmult(xymtxt, 3, samples.size(), zmtx, samples.size(), 1, b);

	double **params = minit(3, 1);
	mmult(ai, 3, 3, b, 3, 1, params);

	mfree(xymtx, samples.size());
	mfree(xymtxt, 3);
	mfree(zmtx, samples.size());
	mfree(a, 3);
	mfree(b, 3);

	for(int r = 0; r < out.rows(); ++r) {
		std::cerr << "Row " << r << " of " << out.rows() << std::endl;
		for(int c = 0; c < out.cols(); ++c) 
			out.set(c, r, params[0][0] + (out.toX(c) - cx) * params[1][0] + (out.toY(r) - cy) * params[2][0]);
	}

	mfree(params, 3);
}

void usage() {
	std::cout << "Usage: rastfit [options]" << std::endl
			<< "This program will produce an adjustment raster to adjust one raster's " << std::endl
			<< "elevations to match another's using a selected algorithm. The mask is " << std::endl
			<< "optional and used to limit the placement of samples to where elevations are " << std::endl
			<< "known to be good. The output is an adjustment raster with the same extent " << std::endl
			<< "as the input raster, with the new resolution. " << std::endl
			<< " -t -- type:        The type of adjustment. " << std::endl
			<< "                    nn     - natural neighbours." << std::endl
			<< "                    planar - plane fit. " << std::endl
			<< " -m -- mask:        an optional raster whose non-zero pixels dictate the locations of allowable samples." << std::endl
			<< " -r -- resolution:  the pixel size of the output." << std::endl
			<< " -s -- samples:     the number of samples." << std::endl
			<< " -b -- base file:   the raster to adjust the adjustment raster to." << std::endl
			<< " -a -- adjustment   file: the raster to adjust." << std::endl
			<< " -o -- output file: the output file." << std::endl;
}

int main(int argc, char **argv) {

 	try {

 		std::string basefile;
 		std::string adjfile;
 		std::string maskfile;
  		std::string outfile;
  		std::string type;
 		unsigned int numSamples = 0;
 		double resolution = 0.0;

 		for(int i = 0; i < argc; ++i) {
 			std::string p(argv[i]);
 			if(p == "-t") {
 				type = argv[++i];
 			} else if(p == "-m") {
 				maskfile = argv[++i];
 			} else if(p == "-s") {
 				numSamples = atoi(argv[++i]);
 			} else if(p == "-o") {
 				outfile = argv[++i];
 			} else if(p == "-a") {
 				adjfile = argv[++i];
 			} else if(p == "-r") {
 				resolution = atof(argv[++i]);
 			} else if(p == "-b") {
 				basefile = argv[++i];
 			}
 		}

 		if(numSamples <= 1)
 			throw "Please try more than one sample.";
 		if(resolution <= 0)
 			throw "Invalid resolution.";
 		if(basefile.empty())
 			throw "Base file is required.";
 		if(outfile.empty())
 			throw "Output file is required.";
 		if(adjfile.empty())
 			throw "Adjustment file is required.";

 		if(type == "nn") {
 			adjustNN(basefile, adjfile, maskfile, outfile, numSamples, resolution);
 		} else if(type == "plane") {
 			adjustPlanar(basefile, adjfile, maskfile, outfile, numSamples, resolution);
 		} else {
 			throw "Unknown type.";
 		}

 		std::cerr << "END" << std::endl;

 	} catch(const char *e) {
 		std::cerr << e << std::endl;
 		usage();
 		return 1;
 	}

 	return 0;
 }



