#include "ogrsf_frmts.h"

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


class ShapeWriter {
private:
	std::string m_filename;
	bool m_needWrite;
	bool m_on = true;
	int m_id = 0;
	std::vector<float> pts;
public:
	ShapeWriter(const char *filename) {
		m_filename = filename;
		m_needWrite = false;
	}
	void on() {
		m_on = true;
	}
	void off() {
		m_on = false;
	}
	void put(Ray_2 &r) {
		if(!m_on) return;
		double x = CGAL::to_double(r.source().x());
		double  y = CGAL::to_double(r.source().y());
		pts.push_back(2);
		pts.push_back(++m_id);
		pts.push_back(x);
		pts.push_back(y);
		pts.push_back(x + CGAL::to_double(r.direction().dx()) * 1000.0);
		pts.push_back(y + CGAL::to_double(r.direction().dy()) * 1000.0);
		m_needWrite = true;
	}
	void put(Segment_2 &r) {
		if(!m_on) return;
		pts.push_back(2);
		pts.push_back(++m_id);
		pts.push_back(CGAL::to_double(r.source().x()));
		pts.push_back(CGAL::to_double(r.source().y()));
		pts.push_back(CGAL::to_double(r.target().x()));
		pts.push_back(CGAL::to_double(r.target().y()));
		m_needWrite = true;
	}
	void put(Segment_2 &r, std::vector<float> &pts) {
		pts.push_back(CGAL::to_double(r.source().x()));
		pts.push_back(CGAL::to_double(r.source().y()));
		pts.push_back(CGAL::to_double(r.target().x()));
		pts.push_back(CGAL::to_double(r.target().y()));
	}
	void put(Point_2 &r) {
		if(!m_on) return;
		pts.push_back(5);
		pts.push_back(++m_id);
		pts.push_back(CGAL::to_double(r.x()));
		pts.push_back(CGAL::to_double(r.y()));
		
		pts.push_back(CGAL::to_double(r.x()+10.0));
		pts.push_back(CGAL::to_double(r.y()));

		pts.push_back(CGAL::to_double(r.x()+10.0));
		pts.push_back(CGAL::to_double(r.y()+10.0));

		pts.push_back(CGAL::to_double(r.x()));
		pts.push_back(CGAL::to_double(r.y()+10.0));

		pts.push_back(CGAL::to_double(r.x()));
		pts.push_back(CGAL::to_double(r.y()));
		m_needWrite = true;
	}
	void put(Iso_rectangle_2 &r) {
		if(!m_on) return;
		float x1 = CGAL::to_double(r.min().x());
		float y1 = CGAL::to_double(r.min().y());
		float x2 = CGAL::to_double(r.max().x());
		float y2 = CGAL::to_double(r.max().y());
		pts.push_back(5);
		pts.push_back(++m_id);
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
		m_needWrite = true;
	}
	void put(Polygon_2 &p) {
		if(!m_on) return;
		Polygon_2::Vertex_const_circulator c = p.vertices_circulator(), done(c);
		std::list<float> tmp;
		Point_2 pt0 = *c;
		tmp.push_back(CGAL::to_double(pt0.x()));
		tmp.push_back(CGAL::to_double(pt0.y()));
		do {
			Point_2 pt = *c;
			tmp.push_back(CGAL::to_double(pt.x()));
			tmp.push_back(CGAL::to_double(pt.y()));
			tmp.push_back(CGAL::to_double(pt.x()));
			tmp.push_back(CGAL::to_double(pt.y()));
		} while(++c != done);
		tmp.push_back(CGAL::to_double(pt0.x()));
		tmp.push_back(CGAL::to_double(pt0.y()));
		pts.push_back((float) tmp.size() / 2);
		pts.push_back(++m_id);
		pts.insert(pts.end(), tmp.begin(), tmp.end());
		m_needWrite = true;
	}
	void write() {
		if(!m_needWrite) return;
		OGRRegisterAll();
		OGRSFDriver *drv = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
		if(!drv)
			throw "Couldn't find shapefile driver.";
		OGRDataSource *ds = drv->CreateDataSource(m_filename.c_str(), NULL);
		if(!ds)
			throw "Failed to open shapefile.";
		OGRLayer *ly = ds->CreateLayer("shapes", NULL, wkbLineString, NULL);
		if(!ly)
			throw "Failed to create layer.";
		OGRFieldDefn fdfn("id", OFTInteger);
		ly->CreateField(&fdfn);
		int count = 0;
		unsigned int i = 0;
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
		m_needWrite = false;
	}
	~ShapeWriter() {
		write();	
	}
};
