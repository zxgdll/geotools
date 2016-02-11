#include "ogrsf_frmts.h"

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
typedef Delaunay::Finite_faces_iterator     DFinite_faces_iterator;
typedef Delaunay::Edge_circulator           DEdge_circulator;
typedef Delaunay::Face_handle               DFace_handle;
typedef Delaunay::Face                      DFace;
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


class ShapeWriter {
private:
	std::string s_filename;
	bool s_needWrite;
	bool s_on = true;
	int s_id = 0;
	std::vector<float> pts;
public:
	ShapeWriter() {
		s_filename = "/tmp/shapes.shp";
		s_needWrite = false;
	}
	void on() {
		s_on = true;
	}
	void off() {
		s_on = false;
	}
	Segment_2 *getHalfedgeAsSegment(VHalfedge &he, Iso_rectangle_2 &bounds) {
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
		CGAL::Object inter = CGAL::intersection(bounds, ray);
		CGAL::Point_2<K> p;
		CGAL::Segment_2<K> s;
		if(CGAL::assign(p, inter)) {
			return new Segment_2(ray.source(), p);
		} else if(CGAL::assign(s, inter)) {
			return new Segment_2(s);			
		} else {
			return NULL;
		}
	}
	void put(Voronoi &v, Iso_rectangle_2 &bounds) {
		++s_id;
		auto it = v.faces_begin();
		do {
			std::vector<float> _pts;
			auto hc = (*it).ccb(), done(hc);
			do {
				if(hc->has_source() && hc->has_target()) {
					Segment_2 seg(hc->source()->point(), hc->target()->point());
					put(seg, _pts);
				} else {
					Segment_2 *seg = getHalfedgeAsSegment(*hc, bounds);
					if(seg != NULL) {
						put(*seg, _pts);
						delete seg;
					}
				}
			} while(++hc != done);

			// Add the boundary vertices that are inside the current face.
			for(int i = 0; i < 4; ++i) {
				VLocate_result lr = v.locate(bounds.vertex(i));
				VFace_handle *fh = boost::get<VFace_handle>(&lr);
				if(fh && **fh == *it) {
					_pts.push_back(bounds.vertex(i).x());
					_pts.push_back(bounds.vertex(i).y());
				}
			}

			pts.push_back(_pts.size() / 2.0);
			pts.push_back(s_id);
			for(auto it2 = _pts.begin(); it2 != _pts.end(); ++it2)
				pts.push_back(*it2);

		} while(++it != v.faces_end());

		s_needWrite = true;
	}
	void put(Delaunay &d) {
		/*
		DFinite_faces_iterator it = d.finite_faces_begin();
		do {
			DFace f = *it;
			DEdge_circulator ec = f.edge_circulator(), ec0 = ec;
			do {
				Segment_2 seg(ec->source().point(), ec->target().point());
				put(seg);
			} while(++ec != ec0);
		} while(it != d.finite_faces_end());
		*/
	}
	void put(Ray_2 &r) {
		if(!s_on) return;
		float x = r.source().x();
		float y = r.source().y();
		pts.push_back(2);
		pts.push_back(++s_id);
		pts.push_back(x);
		pts.push_back(y);
		pts.push_back(x + r.direction().dx() * 1000.0);
		pts.push_back(y + r.direction().dy() * 1000.0);
		s_needWrite = true;
	}
	void put(Segment_2 &r) {
		if(!s_on) return;
		pts.push_back(2);
		pts.push_back(++s_id);
		pts.push_back(r.source().x());
		pts.push_back(r.source().y());
		pts.push_back(r.target().x());
		pts.push_back(r.target().y());
		s_needWrite = true;
	}
	void put(Segment_2 &r, std::vector<float> &pts) {
		pts.push_back(r.source().x());
		pts.push_back(r.source().y());
		pts.push_back(r.target().x());
		pts.push_back(r.target().y());
	}
	void put(Point_2 &r) {
		if(!s_on) return;
		pts.push_back(5);
		pts.push_back(++s_id);
		pts.push_back(r.x());
		pts.push_back(r.y());
		
		pts.push_back(r.x()+10.0);
		pts.push_back(r.y());

		pts.push_back(r.x()+10.0);
		pts.push_back(r.y()+10.0);

		pts.push_back(r.x());
		pts.push_back(r.y()+10.0);

		pts.push_back(r.x());
		pts.push_back(r.y());
		s_needWrite = true;
	}
	void put(Iso_rectangle_2 &r) {
		if(!s_on) return;
		float x1 = r.min().x();
		float y1 = r.min().y();
		float x2 = r.max().x();
		float y2 = r.max().y();
		pts.push_back(5);
		pts.push_back(++s_id);
		pts.push_back(x1);
		pts.push_back(y1);
		
		pts.push_back(x2);
		pts.push_back(y1);
		
		pts.push_back(x2);
		pts.push_back(y2);
		
		pts.push_back(x1);
		pts.push_back(y2);
		
		pts.push_back(x1);
		pts.push_back(y1);
		s_needWrite = true;
	}
	void put(Polygon_2 &p) {
		if(!s_on) return;
		Polygon_2::Vertex_const_circulator c = p.vertices_circulator(), done(c);
		std::list<float> tmp;
		Point_2 pt0 = *c;
		tmp.push_back(pt0.x());
		tmp.push_back(pt0.y());
		do {
			Point_2 pt = *c;
			tmp.push_back(pt.x());
			tmp.push_back(pt.y());
			tmp.push_back(pt.x());
			tmp.push_back(pt.y());
		} while(++c != done);
		tmp.push_back(pt0.x());
		tmp.push_back(pt0.y());
		pts.push_back((float) tmp.size() / 2);
		pts.push_back(++s_id);
		pts.insert(pts.end(), tmp.begin(), tmp.end());
		s_needWrite = true;
	}
	void write() {
		if(!s_needWrite) return;
		OGRRegisterAll();
		OGRSFDriver *drv = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
		if(!drv)
			throw "Couldn't find shapefile driver.";
		OGRDataSource *ds = drv->CreateDataSource(s_filename.c_str(), NULL);
		if(!ds)
			throw "Failed to open shapefile.";
		OGRLayer *ly = ds->CreateLayer("shapes", NULL, wkbLineString, NULL);
		if(!ly)
			throw "Failed to create layer.";
		OGRFieldDefn fdfn("id", OFTInteger);
		ly->CreateField(&fdfn);
		int count = 0;
		int i = 0;
		while(i < pts.size()) {
			OGRFeature *feat = OGRFeature::CreateFeature(ly->GetLayerDefn());
			OGRLineString line;
			count = (int) pts[i++];
			int id = (int) pts[i++];
			feat->SetField(0, id++);
			for(int j = 0; j < count; ++j, i += 2) {
				double x1 = pts[i], y1 = pts[i + 1];
				//std::cerr << x1 << "," << y1 << std::endl;
				line.addPoint(x1, y1, 0.0);
			}
			feat->SetGeometry(&line);
			ly->CreateFeature(feat);
			OGRFeature::DestroyFeature(feat);
		}
		OGRDataSource::DestroyDataSource(ds);
		s_needWrite = false;
	}
	~ShapeWriter() {
		write();	
	}
};