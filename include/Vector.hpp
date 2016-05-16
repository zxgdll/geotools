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

#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Point.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>
 
#include "ogr_spatialref.h"
#include "ogr_geometry.h"
#include "ogrsf_frmts.h"

#include "geotools.h"

using namespace boost::algorithm;

class Vector;

/**
 * Holds information about a newly-added geometry, allowing the
 * caller to add attributes or (potentially) modify the 
 * geometry before it is committed to a vector layer.
 */
class Geom {
friend class Vector;
private:
	OGRFeature *m_feat;
	OGRLayer *m_layer;

	Geom(OGRFeature *feat, OGRLayer *layer) {
		m_feat = feat;
		m_layer = layer;
	}

public:
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
		m_layer->CreateFeature(m_feat);
		OGRFeature::DestroyFeature(m_feat);
	}
};

/**
 * Represents a vector file with a single layer, with a single
 * geometry type that allows the caller to easily create and
 * add geometries with attributes to a data set.
 */
class Vector {
protected:
	OGRDataSource *m_ds;
	OGRLayer *m_layer;
	int m_type;
	std::map<std::string, int> m_atypes;

	OGRwkbGeometryType getGeomType(int type) {
		switch(type) {
		case POINT:
			return wkbPoint;
		case LINE:
			return wkbLineString;
		case POLYGON:
			return wkbPolygon;
		default:
			throw "Unknown or unimplemented geometry type.";
		}
	}
	
public:
	/** Point geometry type. */
	static const int POINT = 1;
	/** Line geometry type. */
	static const int LINE = 2;
	/** MultiLine geometry type. */
	static const int MULTILINE = 3;
	/** Polygon geometry type. */
	static const int POLYGON = 3;
	// TODO: More geometry types.
	
	/** String attribute type. */
	static const int STRING = 1;
	/** Integer attribute type. */
	static const int INTEGER = 2;
	/** Double attribute type. */
	static const int DOUBLE = 3;
	// TODO: More attribute types.
	
	/**
	 * Construct a Vector with the given output file name, geometry 
	 * type and projection information.
	 */
	Vector(const std::string &filename, int type, const std::string &proj, 
		const std::string &vecType = std::string("ESRI Shapefile")) {
		m_type = 0;
		m_layer = nullptr;
		m_ds = nullptr;
		m_type = type;

		OGRRegisterAll();
	
		OGRwkbGeometryType gtype = getGeomType(type);
		
		OGRSpatialReference *gproj = 0;
		if(!proj.empty()) {
			std::string chunk = proj.substr(0, 5);
			to_lower(chunk);
			if(starts_with(chunk, "epsg:")) {
				gproj->importFromEPSG(atoi(proj.substr(5).c_str()));
			} else {
				char *p = (char *) malloc((unsigned long) proj.size() + 1);
				memcpy(p, proj.c_str(), proj.size() + 1);
				gproj = new OGRSpatialReference();
				gproj->importFromWkt(&p);
			}
		}

		OGRSFDriver *drv = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(vecType.c_str());
		if(!drv)
			throw "Failed to load Shapefile driver.";
		if(!(m_ds = drv->CreateDataSource(filename.c_str(), 0)))
			throw "Failed to create vector data source. It may already exist.";
		if(!(m_layer = m_ds->CreateLayer("layer", gproj, gtype, 0)))
			throw "Failed to create vector layer.";
	}

	/**
	 * Construct a Vector with the given output file name, geometry 
	 * type and projection information. Also configures fields from 
	 * the given mapping of attribute properties.
	 */
	Vector(const std::string &filename, int type, const std::string &proj, 
		const std::map<std::string, int> &attributes, const std::string &vecType = std::string("ESRI Shapefile")) :
		Vector(filename, type, proj, vecType) {
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

	/**
	 * Add a point to the geometry with the given coordinates. Returns
	 * a unique_ptr to the Geom object, to which attributes can be added.
	 */
	std::unique_ptr<Geom> addPoint(double x, double y, double z = 0.0) {
		if(m_type != POINT) throw "This is not a point layer.";
		OGRFeature *feat = OGRFeature::CreateFeature(m_layer->GetLayerDefn());
		OGRPoint pt(x, y, z);
		feat->SetGeometry(&pt);
		std::unique_ptr<Geom> ret(new Geom(feat, m_layer));
		return ret;
	}

	std::unique_ptr<Geom> addMultiLine(geos::geom::MultiLineString &line) {
		if(m_type != MULTILINE) throw "This is not a multiline layer.";
		const GEOSContextHandle_t gctx = OGRGeometry::createGEOSContext();
		OGRGeometry *geom = OGRGeometryFactory::createFromGEOS(gctx, (GEOSGeom) &line);
		OGRFeature *feat = OGRFeature::CreateFeature(m_layer->GetLayerDefn());
		feat->SetGeometry(geom);
		std::unique_ptr<Geom> ret(new Geom(feat, m_layer));
		return ret;

	}

	std::unique_ptr<Geom> addMultiLine(const std::vector<std::vector<std::tuple<double, double, double> > > &line) {
		if(m_type != MULTILINE) throw "This is not a multiline layer.";
		OGRFeature *feat = OGRFeature::CreateFeature(m_layer->GetLayerDefn());
		OGRMultiLineString mline;
		
		for(std::vector<std::tuple<double, double, double > > seg:line) {
			OGRLineString lseg;
			int i = 0;
			for(std::tuple<double, double, double> pt:seg)
				lseg.setPoint(i++, std::get<0>(pt), std::get<1>(pt), std::get<2>(pt));
			mline.addGeometry(&lseg);
		}
		
		feat->SetGeometry(&mline);
		std::unique_ptr<Geom> ret(new Geom(feat, m_layer));
		return ret;

	}

	/**
	 * Add a line to the geometry. The argument is a list of tuples 
	 * containing three doubles, x, y, z.
	 * Returns a unique_ptr to the Geom object, to which attributes can be added.
	 */
	std::unique_ptr<Geom> addLine(const std::vector<std::tuple<double, double, double> > &points) {
		if(m_type != LINE) throw "This is not a line layer.";
		OGRFeature *feat = OGRFeature::CreateFeature(m_layer->GetLayerDefn());
		OGRLineString line;
		int i = 0;
		for(std::tuple<double, double, double> pt:points)
			line.setPoint(i++, std::get<0>(pt), std::get<1>(pt), std::get<2>(pt));
		feat->SetGeometry(&line);
		std::unique_ptr<Geom> ret(new Geom(feat, m_layer));
		return ret;
	}
	
	/**
	 * Add a polygon to the geometry. The argument is a list of tuples 
	 * containing three doubles, x, y, z for the exterior ring, and a list of such
	 * lists for holes/islands.
	 * Returns a unique_ptr to the Geom object, to which attributes can be added.
	 */
	std::unique_ptr<Geom> addPolygon(const std::vector<std::tuple<double, double, double> > &extRing, 
																	 std::vector<std::vector<std::tuple<double, double, double > > > &holes) {
		if(m_type != POLYGON) throw "This is not a polygon layer.";
		OGRFeature *feat = OGRFeature::CreateFeature(m_layer->GetLayerDefn());
		OGRPolygon poly;
		
		{
			OGRLinearRing ring;
			int i = 0;
			for(std::tuple<double, double, double> pt:extRing)
				ring.setPoint(i++, std::get<0>(pt), std::get<1>(pt), std::get<2>(pt));
			ring.setPoint(i++, std::get<0>(extRing[0]), std::get<1>(extRing[0]), std::get<2>(extRing[0]));
			poly.addRing(&ring);
		}
		
		for(std::vector<std::tuple<double, double, double > > hole:holes) {
			OGRLinearRing ring;
			int i = 0;
			for(std::tuple<double, double, double> pt:hole)
				ring.setPoint(i++, std::get<0>(pt), std::get<1>(pt), std::get<2>(pt));
			ring.setPoint(i++, std::get<0>(hole[0]), std::get<1>(hole[0]), std::get<2>(hole[0]));
			poly.addRing(&ring);
		}
		
		feat->SetGeometry(&poly);
		std::unique_ptr<Geom> ret(new Geom(feat, m_layer));
		return ret;
	}

	/**
	 * Add a simple polygon to the geometry. The argument is a list of tuples 
	 * containing three doubles, x, y, z for the exterior ring.
	 * Returns a unique_ptr to the Geom object, to which attributes can be added.
	 */
	std::unique_ptr<Geom> addPolygon(const std::vector<std::tuple<double, double, double> > &extRing) {
		std::vector<std::vector<std::tuple<double, double, double > > > holes;
		return addPolygon(extRing, holes);
	}
	
	~Vector() {
		OGRDataSource::DestroyDataSource(m_ds);
	}
};


#endif /* VECTOR_HPP_ */
