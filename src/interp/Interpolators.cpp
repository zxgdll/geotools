#include <iostream>
#include <cmath>

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

#include "QtGui/QApplication"

#include "interp/IDWInterpolator.hpp"
#include "interp/AvgInterpolator.hpp"
#include "interp/PlanarInterpolator.hpp"
#include "interp/NaturalNeighbourInterpolator.hpp"
#include "interp/SimpleKrigingInterpolator.hpp"
#include "interp/krigeplot/KrigePlot.hpp"

namespace interp {

	namespace kriging {

		namespace detail {

		}

		void SimpleKrigingInterpolator::showVariogram(std::list<InterpPoint> &samples) {
			int _argc = argc();
			char **_argv = argv();
			QApplication qa(_argc, _argv);
			QDialog qd;
			interp::kriging::ui::KrigePlot kp;
			kp.setupUi(&qd);
			std::cerr << samples.size() << " samples" << std::endl;
			kp.setSamples(samples);
			qd.show();
			qa.exec();
		}

		void SimpleKrigingInterpolator::interpolate(Raster<float> &out, std::list<InterpPoint > &samples) {
 			showVariogram(samples);
		}
	}

	namespace idw {

		namespace detail {

			/**
			 * Returns the square of the argument.
			 */
			double _sq(double a) {
				return a*a;
			}

			/**
			 * Returns the squared distance between the two points.
			 */
			double _sdist(InterpPoint &a, InterpPoint &b) {
				return _sq(a.x - b.x) + _sq(a.y - b.y);
			}

			/**
			 * Returns the squared distance between the point and the coordinate given by x, y.
			 */
			double _sdist(InterpPoint &a, double x, double y) {
				return _sq(a.x - x) + _sq(a.y - y);
			}

		}

		/**
		 * Performs IDW interpolation. See interp/IDWInterpolator.hpp
		 */
		void IDWInterpolator::interpolate(Raster<float> &out, std::list<InterpPoint > &samples) {
			using namespace detail;
			std::unique_ptr<Block<float>> blk = out.block();
			while(blk->next()) {
				for(int r = blk->startRow(); r < blk->endRow(); ++r) {
					for(int c = blk->startCol(); c < blk->endCol(); ++c) {
						double z = 0.0;
						double t = 0.0;
						for(auto it = samples.begin(); it != samples.end(); ++it)
							t += _sdist(*it, out.toX(c), out.toY(r));
						for(auto it = samples.begin(); it != samples.end(); ++it)
							z += it->diff() * _sdist(*it, out.toX(c), out.toY(r)) / t;
						out.set(c, r, z);
					}
				}
			}
		}

	} // idw

	namespace avg {

		/**
		 * "Interpolates" but doesn't really. Just sets every pixel to the average
		 * of the samples.
		 */
		void AvgInterpolator::interpolate(Raster<float> &out, std::list<InterpPoint > &samples) {
			double z = 0.0;
			for(auto it = samples.begin(); it != samples.end(); ++it)
				z += it->diff();
			z /= samples.size();
			for(int r = 0; r < out.rows(); ++r) {
				for(int c = 0; c < out.cols(); ++c)
					out.set(c, r, z);
			}
		}

	} // avg

	namespace planar {

		namespace detail {

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

		} // detail


		void PlanarInterpolator::interpolate(Raster<float> &out, std::list<InterpPoint > &samples) {
			using namespace detail;

			double **xymtx = minit(samples.size(), 3);
			double **zmtx = minit(samples.size(), 1);
			int i = 0;
			for(auto it = samples.begin(); it != samples.end(); ++it) {
				xymtx[i][0] = 1.0;
				xymtx[i][1] = it->x;
				xymtx[i][2] = it->y;
				zmtx[i][0] = it->diff();
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
				for(int c = 0; c < out.cols(); ++c)
					out.set(c, r, params[0][0] + (out.toX(c) - cx) * params[1][0] + (out.toY(r) - cy) * params[2][0]);
			}

			mfree(params, 3);

		}

	} // planar

	namespace naturalneighbour {

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

		namespace detail {

			double _max(double a, double b) {
				return a < b ? b : a;
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
					double a = atan2((double) CGAL::to_double(dir.dy()), (double) CGAL::to_double(dir.dx()));
					Point_2 pt1 = he.has_source() ? he.source()->point() : he.target()->point();
					Point_2 pt2((double) CGAL::to_double(pt1.x()) + d * cos(a),
							(double) CGAL::to_double(pt1.y()) + d * sin(a));
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

		} // detail

		/**
		 * Performs a natural neighbours interpolation.
		 */
		void NaturalNeighbourInterpolator::interpolate(Raster<float> &out, std::list<InterpPoint > &samples) {
			using namespace detail;

			// Build a boundary for clipping the voronoi.
			std::list<Point_2> pts;
			pts.push_back(Point_2(out.maxx() + 1000.0, out.maxy() + 1000.0));
			pts.push_back(Point_2(out.maxx() + 1000.0, out.miny() - 1000.0));
			pts.push_back(Point_2(out.minx() - 1000.0, out.miny() - 1000.0));
			pts.push_back(Point_2(out.minx() - 1000.0, out.maxy() + 1000.0));
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
				diffs[id] = pt->diff();
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
				//std::cerr << "Row " << r << " of " << out.rows() << std::endl;
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

	} // naturalneighbour

} // interp


