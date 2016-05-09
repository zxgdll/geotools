/*
 * Vector.hpp
 *
 *  Created on: May 7, 2016
 *      Author: rob
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include <map>
#include <memory>

#include "gdal/ogr_spatialref.h"
#include "gdal/ogr_geometry.h"
#include "gdal/ogrsf_frmts.h"

class Geom {
private:
	OGRFeature *m_feat;
public:
	Geom(OGRFeature *feat) {
		m_feat = feat;
	}

	void setAttribute(const char *name, int value) {
		m_feat->SetField(name, value);
	}

	void setAttribute(const char *name, double value) {
		m_feat->SetField(name, value);
	}

	void setAttribute(const char *name, std::string value) {
		m_feat->SetField(name, value.c_str());
	}

	~Geom() {
		OGRFeature::DestroyFeature(m_feat);
	}
};

class Vector {
protected:
	OGRDataSource *m_ds;
	OGRLayer *m_layer;
	int m_type;
	std::map<std::string, int> m_atypes;

public:
	static const int POINT = 1;
	static const int LINE = 2;
	static const int POLYGON = 3;

	static const int STRING = 1;
	static const int INTEGER = 2;
	static const int DOUBLE = 3;

	Vector(std::string &filename, int type, std::string &proj) {
		m_type = 0;
		m_layer = nullptr;
		m_ds = nullptr;

		OGRRegisterAll();

		OGRwkbGeometryType gtype = wkbNone;
		switch(type) {
		case POINT:
			gtype = wkbPoint;
			break;
		case LINE:
			gtype = wkbLineString;
			break;
		case POLYGON:
			gtype = wkbPolygon;
			break;
		default:
			throw "Unknown geometry type.";
		}

		m_type = type;

		OGRSpatialReference *gproj = 0;
		if(!proj.empty()) {
			char *p = (char *) malloc((unsigned long) proj.size() + 1);
			memcpy(p, proj.c_str(), proj.size() + 1);
			gproj = new OGRSpatialReference();
			gproj->importFromWkt(&p);
		}

		OGRSFDriver *drv = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
		if(!drv)
			throw "Failed to load Shapefile driver.";
		if(!(m_ds = drv->CreateDataSource(filename.c_str(), 0)))
			throw "Failed to create vector data source. It may already exist.";
		if(!(m_layer = m_ds->CreateLayer("layer", gproj, gtype, 0)))
			throw "Failed to create vector layer.";

		OGRFieldDefn idField("gid", OFTInteger);
		m_layer->CreateField(&idField);
	}

	Vector(std::string &filename, int type, std::string &proj, std::map<std::string, int> &attributes) :
		Vector(filename, type, proj) {
		for(auto it = attributes.begin(); it != attributes.end(); ++it) {
			OGRFieldType type;
			m_atypes[it->first] = it->second;
			switch(it->second) {
			case STRING:
				type = OFTString;
				break;
			case INTEGER:
				type = OFTInteger;
				break;
			case DOUBLE:
				type = OFTReal;
				break;
			default:
				throw "Unknown field type.";
			}
			OGRFieldDefn field(it->first.c_str(), type);
			m_layer->CreateField(&field);
		}
	}

	std::unique_ptr<Geom> addPoint(double x, double y, double z = 0.0) {
		if(m_type != POINT) throw "This is not a point layer.";
		OGRFeature *feat = OGRFeature::CreateFeature(m_layer->GetLayerDefn());
		OGRPoint pt(x, y, z);
		feat->SetGeometry(&pt);
		m_layer->CreateFeature(feat);
		std::unique_ptr<Geom> ret(new Geom(feat));
		return ret;
	}

	~Vector() {
		OGRDataSource::DestroyDataSource(m_ds);
	}
};


#endif /* VECTOR_HPP_ */
