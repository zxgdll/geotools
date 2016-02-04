#include "ogrsf_frmts.h"
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<float> 							K;
typedef CGAL::Iso_rectangle_2<K>														Iso_rectangle_2;
typedef CGAL::Polygon_2<K>																Polygon_2;
typedef K::Point_2 							Point_2;
typedef K::Ray_2							Ray_2;
typedef K::Segment_2						Segment_2;
typedef K::Line_2							Line_2;

class ShapeWriter {
private:
	std::string s_filename;
	bool s_needWrite;
	bool s_on = true;
public:
	std::vector<float> pts;
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
	void put(Ray_2 &r) {
		if(!s_on) return;
		float x = r.source().x();
		float y = r.source().y();
		pts.push_back(2);
		pts.push_back(x);
		pts.push_back(y);
		pts.push_back(x + r.direction().dx() * 1000.0);
		pts.push_back(y + r.direction().dy() * 1000.0);
		s_needWrite = true;
	}
	void put(Segment_2 &r) {
		if(!s_on) return;
		pts.push_back(2);
		pts.push_back(r.source().x());
		pts.push_back(r.source().y());
		pts.push_back(r.target().x());
		pts.push_back(r.target().y());
		s_needWrite = true;
	}
	void put(Point_2 &r) {
		if(!s_on) return;
		pts.push_back(5);
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
		int count = 0;
		int i = 0;
		while(i < pts.size()) {
			OGRFeature *feat = OGRFeature::CreateFeature(ly->GetLayerDefn());
			OGRLineString line;
			count = (int) pts[i++];
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