#include <iostream>
#include <cmath>
#include <vector>
#include <list>

#include <omp.h> 

#ifdef WITH_CGAL
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
#endif

#include "nanoflann.hpp"

#ifdef WITH_QT
#include "QtGui/QApplication"

#include <interp/ui/KrigePlot.hpp>
#endif

#include <Eigen/LU>

#include "interp/IDWInterpolator.hpp"
#include "interp/AvgInterpolator.hpp"
#include "interp/PlanarInterpolator.hpp"
#include "interp/NaturalNeighbourInterpolator.hpp"
#include "interp/SimpleKrigingInterpolator.hpp"

namespace interp {

	namespace detail {

		// See nanoflann examples: https://github.com/jlblancoc/nanoflann/blob/master/examples/pointcloud_kdd_radius.cpp#L119
		struct PointCloud
		{
			std::vector<InterpPoint> pts;

			// Must return the number of data points
			inline size_t kdtree_get_point_count() const { return pts.size(); }

			// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
			inline double kdtree_distance(const double *p1, const size_t idx_p2,size_t /*size*/) const
			{
				const double d0=p1[0]-pts[idx_p2].x;
				const double d1=p1[1]-pts[idx_p2].y;
				//const double d2=p1[2]-pts[idx_p2].z;
				return d0*d0+d1*d1;//+d2*d2;
			}

			// Returns the dim'th component of the idx'th point in the class:
			// Since this is inlined and the "dim" argument is typically an immediate value, the
			//  "if/else's" are actually solved at compile time.
			inline double kdtree_get_pt(const size_t idx, int dim) const
			{
				if (dim==0) return pts[idx].x;
				else if (dim==1) return pts[idx].y;
				else return pts[idx].z;
			}

			// Optional bounding-box computation: return false to default to a standard bbox computation loop.
			//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
			//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
			template <class BBOX>
			bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

		};

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

	namespace kriging {

		void SimpleKrigingInterpolator::computeVariogram(std::list<InterpPoint> &samples,
				std::list<VariogramPoint> &variogram) {
			using namespace interp::detail;

			for(auto it0 = samples.begin(); it0 != samples.end(); ++it0) {
				for(auto it1 = samples.begin(); it1 != samples.end(); ++it1) {
					if(it0->equals(*it1)) continue;
					double dist = sqrt(_sq(it0->x - it1->x) + _sq(it0->y - it1->y));
					if(std::isnan(dist)) {
						std::cerr << it0->x << ", " << it1->x << ":" << it0->y << ", " << it1->y << std::endl;
						std::cerr << _sq(it0->x - it1->x) << ":" << _sq(it0->y - it1->y) << std::endl;
						throw "Bad distance.";
					}
					double diff = _sq(it0->z - it1->z) / 2.0;
					variogram.push_back(VariogramPoint(dist, diff));
				}
			}
		}

		void SimpleKrigingInterpolator::showVariogram(interp::kriging::KrigeArgs &kargs, std::list<VariogramPoint> &variogram) {
#ifdef WITH_QT
			int _argc = argc();
			char **_argv = argv();
			QApplication qa(_argc, _argv);
			QDialog qd;
			interp::ui::KrigePlot kp(&kargs);
			kp.setupUi(&qd);
			std::cerr << variogram.size() << " variogram points" << std::endl;
			kp.setVariogram(variogram);
			qd.show();
			qa.exec();
#else
			throw "Must be compiled with Qt for Kriging.";
#endif
		}

		void SimpleKrigingInterpolator::interpolate(Raster<float> &out, std::list<InterpPoint > &samples) {
			using namespace interp::detail;
			using namespace Eigen;

			std::list<VariogramPoint> variogram;
			interp::kriging::KrigeArgs kargs;
			computeVariogram(samples, variogram);
 			showVariogram(kargs, variogram);
 			if(kargs.status == 1) {
 				std::cerr << "Cancelled" << std::endl;
 				return;
 			}

 			std::cerr << "Params: nugget: " << kargs.nugget << std::endl
					<< "range: " << kargs.range << std::endl
					<< "sill: " << kargs.sill << std::endl;

 			// All-pairs distance matrix.
 			Matrix<double, Dynamic, Dynamic> D(samples.size(), samples.size());
 			// Semivariance matrix.
 			Matrix<double, Dynamic, Dynamic> A(samples.size() + 1, samples.size() + 1);
 			int r = 0, c = 0;
 			for(auto it0 = samples.begin(); it0 != samples.end(); ++it0) {
 				c = 0;
 				for(auto it1 = samples.begin(); it1 != samples.end(); ++it1) {
 					double dist = sqrt(_sq(it0->x - it1->x) + _sq(it0->y - it1->y));
 					D(c, r) = dist;
 					// TODO: Use functor for this?
 					A(c, r) = kargs.model(dist, kargs.nugget, kargs.sill, kargs.range);
 					++c;
 				}
 				++r;
 			}
 			// Pad out the matrix with 1s and a 0.
 			for(unsigned int i = 0; i < samples.size(); ++i) {
 				A(i, A.rows()) = 1;
 				A(A.cols(), i) = 1;
 			}
 			A(A.cols(), A.rows()) = 0;

 			// TODO: Use cell centroid rather than corner for distance in all methods.
			// Grid for distances from cell to control points.
			Matrix<double, Dynamic, Dynamic> b(samples.size(), 1);
			// Grid for modeled variances from cell to control points.
			Matrix<double, Dynamic, Dynamic> d(samples.size() + 1, 1);
 			for(int r = 0; r < out.rows(); ++r) {
 				for(int c = 0; c < out.cols(); ++c) {
 					int i = 0;
 					for(auto it = samples.begin(); it != samples.end(); ++it) {
 						double dist = sqrt(_sq(out.toX(c) - it->x) + _sq(out.toY(r) - it->y));
 						d(i, 0) = dist;
 						b(i, 0) = kargs.model(dist, kargs.nugget, kargs.sill, kargs.range);
 						++i;
 					}
 					d(samples.size(), 0) = 1.0; // Last row is 1
 					Matrix<double, Dynamic, Dynamic> Ai = A.inverse();
 					Matrix<double, Dynamic, Dynamic> w = Ai * b; // Last row is the LeGrangian: ignore.
 					double px = 0.0;
 					i = 0;
 					for(auto it = samples.begin(); it != samples.end(); ++it) {
 						px += it->z * w(i, 0);
 						++i;
 					}
 					out.set(c, r, (float) px);
 				}
 			}
		}
	}

	namespace idw {

		/**
		 * Performs IDW interpolation. See interp/IDWInterpolator.hpp
		 */
		void IDWInterpolator::interpolate(Raster<float> &out, std::list<InterpPoint > &samples) {
			using namespace interp::detail;
			using namespace detail;
			using namespace nanoflann;

			// If we're using all points.
			if(m_neighbours == 0) {

				Block<float> blk = out.block();

				#pragma omp parallel
				{

				// Use a block instance for writing output in a thread-safe way.
				Grid<float> grid;
				int sr = 0, sc = 0;

				while(true) {

					bool run = true;
					#pragma omp critical
					{
						// Try to advance the block. If it fails, quit the loop.
						// Otherwise, get the limits.
						run = blk.next();
					} // omp

					if(run) {
						sr = blk.startRow(), sc = blk.startCol();
						int rows = blk.rows(), cols = blk.cols();
						grid.init(cols, rows);

						// Iterate over the cells in the 
						for(int r = 0; r < rows; ++r) {
							std::cerr << "row " << (r + sr) << std::endl;
							for(int c = 0; c < cols; ++c) {
								double z = 0.0;
								double t = 0.0;
								for(auto it = samples.begin(); it != samples.end(); ++it) {
									double d = pow(_sdist(*it, out.toX(c), out.toY(r)), m_exponent);
									if(d < std::numeric_limits<double>::lowest()) {
										z = it->z;
										t = 1.0;
										break;
									} else {
										z += it->z / d;
										t += 1 / d;
									}
								}
								grid(c, r, z / t);
							}
						}
					}

					#pragma omp critical
					{
						out.set(sc, sr, grid);
					} // omp

				}
				} // omp

			// If we're using a subset of neighbours.
			// TODO: NN, or radius?
			} else {

				// Set the number of neighbours to use. If the config is for <= 0, use them all.
				// Otherwise, pick the smaller of sample size and neighbours.
				const int num = m_neighbours <= 0 ? samples.size() : (samples.size() >= m_neighbours ? m_neighbours : samples.size());

				// Prepare a kdtree to find neighbours for computing idw.
				PointCloud pc;
				pc.pts.assign(samples.begin(), samples.end());
				typedef KDTreeSingleIndexAdaptor<
					L2_Simple_Adaptor<double, PointCloud>,
					PointCloud,
					2
				> kd_tree;
				kd_tree index(2, pc, KDTreeSingleIndexAdaptorParams(10)); // TODO: Best leaf size?
				index.buildIndex();

				std::vector<unsigned long> idx(num); 	// Point indices
				std::vector<double> dist(num);			// Distance from query point.

				// Get the block to navigate the raster's internal block structure.
				Block<float> blk = out.block();

				#pragma omp parallel
				{

				// Use a block instance for writing output in a thread-safe way.
				Grid<float> grid;
				int sr = 0, sc = 0;

				while(true) {

					bool run = true;
					#pragma omp critical
					{
						// Try to advance the block. If it fails, quit the loop.
						// Otherwise, get the limits.
						run = blk.next();
					} // omp

					if(run) {
						sr = blk.startRow(), sc = blk.startCol();
						int rows = blk.rows(), cols = blk.cols();
						grid.init(cols, rows);

						// Iterate over the cells in the
						for(int r = 0; r < rows; ++r) {
							std::cerr << "row " << (r + sr) << std::endl;
							for(int c = 0; c < cols; ++c) {
								// Create a query point.
								const double query[2] = {out.toX(c + sc), out.toY(r + sr)};
								// Find the results in the knn tree.
								index.knnSearch(query, num, &idx[0], &dist[0]);
								double z = 0.0;
								double t = 0.0;
								for(int i = 0; i < num; ++i) {
									InterpPoint pt = pc.pts[idx[i]];
									double d = pow(dist[i], m_exponent);
									if(d <= std::numeric_limits<float>::lowest()) {
										// If distance is close to zero, use the value of the
										// nearest point.
										z = pt.z;
										t = 1;
										break;
									} else {
										z += pt.z / d;
										t += 1 / d;
									}
								}
								grid(c, r, z / t);
							}
						}

						#pragma omp critical
						{
							out.set(sc, sr, grid);
						} // omp

					}
				}

				} // omp
			}
		}

	} // idw

	namespace avg {

		/**
		 * "Interpolates" but doesn't really. Just sets every pixel to the average
		 * of the samples.
		 */
		void AvgInterpolator::interpolate(Raster<float> &out, std::list<InterpPoint > &samples) {

			if(m_neighbours <= 0 || m_neighbours >= samples.size()) {

				double z = 0.0;
				for(auto it = samples.begin(); it != samples.end(); ++it)
					z += it->z;
				z /= samples.size();
				for(int r = 0; r < out.rows(); ++r) {
					for(int c = 0; c < out.cols(); ++c)
						out.set(c, r, z);
				}

			} else {

				using namespace interp::detail;
				using namespace detail;
				using namespace nanoflann;

				// Prepare a kdtree to find neighbours for computing idw.
				PointCloud pc;
				pc.pts.assign(samples.begin(), samples.end());
				typedef KDTreeSingleIndexAdaptor<
					L2_Simple_Adaptor<double, PointCloud>,
					PointCloud,
					2
				> kd_tree;
				kd_tree index(2, pc, KDTreeSingleIndexAdaptorParams(10)); // TODO: Best leaf size?
				index.buildIndex();

				std::vector<unsigned long> idx(m_neighbours); 	// Point indices
				std::vector<double> dist(m_neighbours);			// Distance from query point.

				Block<float> blk = out.block();
				while(blk.next()) {
					for(int r = blk.startRow(); r < blk.endRow(); ++r) {
						std::cerr << "row " << r << std::endl;
						for(int c = blk.startCol(); c < blk.endCol(); ++c) {
							const double query[2] = {out.toX(c), out.toY(r)};
							index.knnSearch(query, m_neighbours, &idx[0], &dist[0]);
							double z = 0.0;
							for(unsigned int i = 0; i < m_neighbours; ++i)
								z += pc.pts[idx[i]].z;
							out.set(c, r, z / m_neighbours);
						}
					}
				}
			}
		}

	} // avg

	namespace planar {

		namespace detail {

			/**
			 * Find the centroid of a point cloud.
			 */
			void mcentroid(Eigen::Matrix<double, Eigen::Dynamic, 3> &mtx,
					int rows, int xcol, int ycol, double *x, double *y) {
				*x = 0, *y = 0;
				for(int r = 0; r < rows; ++r) {
					*x += mtx(r, xcol);
					*y += mtx(r, ycol);
				}
				*x /= rows;
				*y /= rows;
				for(int r = 0; r < rows; ++r) {
					mtx(r, xcol) -= *x;
					mtx(r, ycol) -= *y;
				}
			}

			/**
			 * Print a matrix.
			 */
			/*
			void mprint(double **mtx, int rows, int cols) {
				std::cerr << "mtx " << rows << "," << cols << std::endl;
				for(int r = 0; r < rows; ++r) {
					for(int c = 0; c < cols; ++c) {
						std::cerr << mtx[r][c] << " ";
					}
					std::cerr << std::endl;
				}
			}
			*/

			/**
			 * Compute the planar function parameters.
			 * params is the parameter list
			 * cx and cy are the centroid coordinates; use these to offset the arguments.
			 */
			void computeParams(std::list<interp::InterpPoint> &samples,
					Eigen::Matrix<double, 3, 1> &params, double *cx, double *cy) {
				using namespace Eigen;
				// x,y matrix
				Matrix<double, Dynamic, 3> xymtx(samples.size(), 3);
				// z matrix
				Matrix<double, Dynamic, 1> zmtx(samples.size(), 1);
				int i = 0;
				for(auto it = samples.begin(); it != samples.end(); ++it) {
					xymtx(i, 0) = 1.0;
					xymtx(i, 1) = it->x;
					xymtx(i, 2) = it->y;
					zmtx(i, 0) = it->z;
					++i;
				}

				// find the centroid of the xy matrix -- this supplies
				// an offset to help produce a non-zero determinant.
				mcentroid(xymtx, samples.size(), 1, 2, cx, cy);

				Matrix<double, 3, Dynamic> xymtxt = xymtx.transpose();
				Matrix<double, 3, 3> a = xymtxt * xymtx;
				Matrix<double, 3, 3> ai = a.inverse();
				Matrix<double, Dynamic, Dynamic> b = xymtxt * zmtx;
				Matrix<double, Dynamic, Dynamic> _params = ai * b;
				for(int i = 0; i < 3; ++i)
					params(i, 0) = _params(i, 0);
			}

			/**
			 * This utility method interpolates for a single coordinate, rather than
			 * a raster, and returns the interpolated value.
			 */
			double interpolate(double x, double y, std::list<interp::InterpPoint > &samples) {
				using namespace detail;
				using namespace Eigen;
				Matrix<double, 3, 1> params(3, 1);
				double cx, cy; // Centroid for offset.
				computeParams(samples, params, &cx, &cy);
				double z = params(0) + (x - cx) * params(1) + (y - cy) * params(1);
				return z;
			}

		} // detail


		void PlanarInterpolator::interpolate(Raster<float> &out, std::list<interp::InterpPoint > &samples) {
			using namespace detail;
			using namespace Eigen;
			Matrix<double, 3,1 > params(3, 1);
			double cx, cy; // Centroid for offset.
			computeParams(samples, params, &cx, &cy);
			for(int r = 0; r < out.rows(); ++r) {
				for(int c = 0; c < out.cols(); ++c)
					out.set(c, r, (float) (params(0) + (out.toX(c) - cx) * params(1) + (out.toY(r) - cy) * params(2)));
			}
		}

	} // planar

	namespace naturalneighbour {

#ifdef WITH_CGAL
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
#endif
		/**
		 * Performs a natural neighbours interpolation.
		 */
		void NaturalNeighbourInterpolator::interpolate(Raster<float> &out, std::list<InterpPoint > &samples) {

#ifdef WITH_CGAL
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
				diffs[id] = pt->z;
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
#else
			throw "Must be compiled with CGAL for natural neighbours.";
#endif
		}

	} // naturalneighbour

} // interp


