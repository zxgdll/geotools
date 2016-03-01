/*
 * lasstats computes zonal statistics on a point
 * cloud for polygons defined in a Shapefile, or grid cells in a raster.
 *
 * If a shapefile is given, each polygon is updated with attributes for mean,
 * max, min, median, stddev and variance and, optionally, quantiles.
 *
 * If a raster is given, a table is generated giving the statistics for each
 * unique cell ID and point class. That is, if the raster has cells with ID 1,
 * 2 and 3, and the point cloud has points classified 1 and 2, the table will
 * have 6 rows.
 *
 *  Author: Rob Skelly rob@dijital.ca
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <list>
#include <climits>
#include <memory>
#include <iomanip>
#include <cmath>

#include <ogr_api.h>
#include <ogrsf_frmts.h>
#include <gdal_priv.h>

#include <boost/algorithm/string/predicate.hpp>

#include <geos/geom/GeometryCollection.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/Point.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CoordinateSequenceFactory.h>

#include <liblas/liblas.hpp>
 
#include "Util.hpp"

namespace las = liblas;
namespace geom = geos::geom;
namespace alg = boost::algorithm;

#define RASTER 1
#define VECTOR 2
#define CSV 3

/**
 * Squares a double.
 */
double _sq(double a) {
	return a * a;
}

/**
 * Represents a single polygon in the input shape file.
 * An instance of this class can comute all the necessary
 * statistics for its polygon.
 */
class Stat {

private:
	// A pointer to the geometry.
	geom::Geometry *m_geom;
	// A reference to the feature containing the geometry and attributes.
	OGRFeature *m_feat;
	// The list of point z values.
	std::vector<double> m_values;
	// Point class.
	int m_ptClass;
	// Cell ID (for rasters).
	int m_cellId;

	bool m_sorted;
	double m_sum;
	double m_min;
	double m_max;
	double m_mean;
	double m_median;
	double m_stddev;
	double m_variance;
	bool m_reset;

	void reset() {
		if(m_reset) {
			m_sorted = false;
			m_sum = nan("");
			m_min = nan("");
			m_max = nan("");
			m_mean = nan("");
			m_median = nan("");
			m_stddev = nan("");
			m_variance = nan("");
			m_reset = false;
		}
	}

	void init() {
		m_feat = nullptr;
		m_geom = nullptr;
		m_ptClass = 0;
		m_cellId = 0;
		m_reset = true;
		reset();
	}
public:

	/**
	 * Construct a Stat object for the point class
	 * using the given geometry and feature.
	 */
	Stat(int ptClass, geom::Geometry *geom, OGRFeature *feat) {
		init();
		m_ptClass = ptClass;
		m_geom = geom; //->clone();
		m_feat = feat; //->Clone();
	}

	/**
	 * Construct a Stat object using the given point class
	 * and cell ID.
	 */
	Stat(int ptClass, int cellId) {
		init();
		m_ptClass = ptClass;
		m_cellId = cellId;
	}

	static int numBands(int numQuantiles) {
		return 8 + numQuantiles + 2; // plus 2 for the 0th and nth.
	}

	static void bandNames(std::vector<std::string> &names, int numQuantiles) {
		std::string n[8] = {"count", "sum", "min", "max", "mean", "median", "variance", "stddev"};
		for(int i = 0; i < 8; ++i)
			names.push_back(n[i]);
		for(int i = 0; i < numQuantiles + 2; ++i)
			names.push_back("q" + std::to_string(i));
	}

	void bandValues(std::vector<double> &values, int numQuantiles) {
		values.push_back(this->count());
		values.push_back(this->sum());
		values.push_back(this->min());
		values.push_back(this->max());
		values.push_back(this->mean());
		values.push_back(this->median());
		values.push_back(this->variance());
		values.push_back(this->stddev());
		double q[numQuantiles];
		quantiles(q, numQuantiles);
		for(int i = 0; i < numQuantiles + 2; ++i)
			values.push_back(q[i]);
	}

	/**
	 * Sorts the values.
	 */
	void sort() {
		if(!m_sorted) {
			std::sort(m_values.begin(), m_values.end());
			m_sorted = true;
		}
	}

	/**
	 * Add a point value to the list.
	 */
	void add(double value) {
		m_values.push_back(value);
		reset();
	}

	/**
	 * Return the number of values.
	 */
	int count() {
		return m_values.size();
	}

	double min() {
		if(std::isnan(m_min) && count() > 0) {
			m_reset = true;
			sort();
			m_min = m_values[0];
		}
		return m_min;
	}

	double max() {
		if(std::isnan(m_max) && count() > 0) {
			m_reset = true;
			sort();
			m_max = m_values[m_values.size() - 1];
		}
		return m_max;
	}

	/**
	 * Return the sum of values.
	 */
	double sum() {
		if(std::isnan(m_sum) && count() > 0) {
			m_reset = true;
			m_sum = 0.0;
			for(auto it = m_values.begin(); it != m_values.end(); ++it)
				m_sum += (*it);
		}
		return m_sum;
	}

	/**
	 * Return the mean of values.
	 */
	double mean () {
		if(std::isnan(m_mean)) {
			m_reset = true;
			m_mean = sum() / count();
		}
		return m_mean;
	}

	/**
	 * Return the median of values.
	 */
	double median() {
		int num;
		if(std::isnan(m_median) && (num = count()) > 0) {
			m_reset = true;
			if(num > 0) {
				sort();
				if(num % 2 == 0) {
					m_median = (m_values[num / 2] + m_values[num / 2 - 1]) / 2.0;
				} else {
					m_median = m_values[num / 2];
				}
			}
		}
		return m_median;
	}

	/**
	 * Return the variance of values.
	 */
	double variance() {
		if(std::isnan(m_variance)) {
			m_reset = true;
			double m = mean();
			double s = 0.0;
			for(auto it = m_values.begin(); it != m_values.end(); ++it)
				s += _sq((*it) - m);
			m_variance = s / count();
		}
		return m_variance;
	}

	/**
	 * Return the standard deviation of values.
	 */
	double stddev() {
		if(std::isnan(m_stddev)) {
			m_reset = true;
			m_stddev = sqrt(variance());
		}
		return m_stddev;
	}

	/**
	 * Return the given number of quantiles.
	 * This populates the given array with n+1
	 * values, so a 4-quantile has five values corresponding
	 * to the lower and upper limits, and the three
	 * quartiles.
	 */
	void quantiles(double *q, int num) {
		if(num < 1)
			throw "Less than one quantile doesn't make sense.";
		int cnt = count();
		m_reset = true;
		for(int i = 0; i < num + 2; ++i) {
			if(cnt == 0) {
				q[i] = nan("");
			} else {
				int c = (int) ceil(((double) i / (num + 1)) * (cnt - 1));
				q[i] = m_values[c];
			}
		}
	}

	geom::Geometry *geom() {
		return m_geom;
	}

	OGRFeature *feat() {
		return m_feat;
	}

	~Stat() {
		//if(m_geom != nullptr)
		//	delete m_geom;
		//if(m_feat != nullptr)
		//	OGRFeature::DestroyFeature(m_feat);
	}


};

void doVector(std::string &shapefile, std::string &outfile, std::string &layername,
		std::vector<std::string> &files, std::vector<int> &classes, int numQuantiles, int outType) {

	const GEOSContextHandle_t gctx = OGRGeometry::createGEOSContext();
	const geom::GeometryFactory *gf = geom::GeometryFactory::getDefaultInstance();
	const geom::CoordinateSequenceFactory *cf = gf->getCoordinateSequenceFactory();

	OGRLayer *layer;
	OGRFeature *feat;
	OGRFeatureDefn *featDefn;
	OGRwkbGeometryType type;
	geom::Geometry *geom;
	geom::Geometry *roi = NULL;
	std::vector<Stat> stats;

	std::cerr << "Opening " << shapefile << std::endl;

 	OGRRegisterAll();

	// Load the shapefile and get the polygons from it.
	OGRDataSource *srcDs = OGRSFDriverRegistrar::Open(shapefile.c_str(), FALSE);
	if(srcDs == NULL)
		throw "Couldn't open shapefile.";

	if(layername.empty()) {
		layer = srcDs->GetLayer(0);
	} else {
		layer = srcDs->GetLayerByName(layername.c_str());
	}
	if(layer == NULL)
		throw "Couldn't get layer.";

	type = layer->GetGeomType();
	if(type != wkbPolygon)
		throw "Geometry must be polygon.";

	featDefn = layer->GetLayerDefn();

	// For each polygon, create a Stat instance.
	// Also create a unified geometry representing the bounds of
	// all polygons.
	while((feat = layer->GetNextFeature()) != NULL) {
		// Export the geometry for querying.
		geom = (geom::Geometry *) feat->GetGeometryRef()->exportToGEOS(gctx);
		// Create a unioned geom to filter las files.
		if(roi == NULL) {
			roi = geom;
		} else {
			roi = roi->Union(geom);
		}
		// Add a new Stat object.
		// TODO: Point classes for vectors. Setting to 0 for now.
		stats.push_back(Stat(0, geom, feat));
	}

	// If there are no geometries, just quit.
	if(stats.size() == 0)
		throw "No geometries were found.";

	// Iterate over LAS files, process each one.
	las::ReaderFactory rf;
	std::cerr << "Processing files." << std::endl;
	for(unsigned int i = 0; i < files.size(); ++i) {

		std::ifstream in(files[i].c_str(), std::ios::in | std::ios::binary);
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();

		std::cout << "Processing file " << files[i] << std::endl;

		// Build the boundary polygon.
		std::vector<geom::Coordinate> coords;
		coords.push_back(geom::Coordinate(h.GetMinX(), h.GetMinY()));
		coords.push_back(geom::Coordinate(h.GetMaxX(), h.GetMinY()));
		coords.push_back(geom::Coordinate(h.GetMaxX(), h.GetMaxY()));
		coords.push_back(geom::Coordinate(h.GetMinX(), h.GetMaxY()));
		coords.push_back(geom::Coordinate(h.GetMinX(), h.GetMinY()));
		geom::CoordinateSequence *cs = cf->create(&coords);
		geom::LinearRing *lr = gf->createLinearRing(cs);
		geom::Polygon *bounds = gf->createPolygon(lr, NULL);

		// If the boundary intersectes the region of interest..
		if(roi->intersects(bounds)) {
			std::cerr << "...file intersects." << std::endl;
			while(r.ReadNextPoint()) {
				las::Point pt = r.GetPoint();
				// ... and this point has a class of interest.
				if(classes.size() != 0 && !Util::inList(classes, pt.GetClassification().GetClass()))
					continue;
				const geom::Coordinate c(pt.GetX(), pt.GetY(), pt.GetZ());
				geom::Point *p = gf->createPoint(c);
				// Add it to each Stat whose polygon contains it.
				for(unsigned int j = 0; j < stats.size(); ++j) {
					if(stats[j].geom()->contains(p)) {
						stats[j].add(c.z);
					}
				}
				delete p;
			}
		}
		in.close();
	}

	// Start building the output shapefile.
	std::cerr << "Updating shapefile." << std::endl;
	OGRRegisterAll();
	OGRSFDriver *drv = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	if(drv == NULL)
		throw "Shapefile driver is not available.";

	OGRDataSource *dstDs = drv->CreateDataSource(outfile.c_str(), NULL);
	if(dstDs == NULL)
		throw "Couldn't create shapefile.";

	layer = dstDs->CreateLayer("stats", NULL, wkbPolygon, NULL);
	if(layer == NULL)
		throw "Couldn't create layer.";

	// Copy the field definitions from the original layer.
	for(int i = 0; i < featDefn->GetFieldCount(); ++i) {
		if(layer->CreateField(featDefn->GetFieldDefn(i)) != OGRERR_NONE)
			throw "Failed to create field."; // + featDefn->GetFieldDefn(i)->GetNameRef();
	}

	// Copy new stats fields to the shapefile.
	OGRFieldDefn fstddev("stddev", OFTReal);
	if(layer->CreateField(&fstddev) != OGRERR_NONE)
		throw "Failed to create stddev field.";

	OGRFieldDefn fmean("mean", OFTReal);
	if(layer->CreateField(&fmean) != OGRERR_NONE)
		throw "Failed to create mean field.";

	OGRFieldDefn fmedian("median", OFTReal);
	if(layer->CreateField(&fmedian) != OGRERR_NONE)
		throw "Failed to create median field.";

	OGRFieldDefn fcount("count", OFTInteger);
	if(layer->CreateField(&fcount) != OGRERR_NONE)
		throw "Failed to create count field.";

	OGRFieldDefn fsum("sum", OFTReal);
	if(layer->CreateField(&fsum) != OGRERR_NONE)
		throw "Failed to create sum field.";

	if(numQuantiles > 0) {
		for(int qt = 0; qt < numQuantiles + 2; ++qt) {
			OGRFieldDefn fqt(("qt" + std::to_string(qt)).c_str(), OFTReal);
			if(layer->CreateField(&fqt) != OGRERR_NONE)
				throw "Failed to create quantile field.";
		}
	}

	double quantiles[numQuantiles > 0 ? numQuantiles + 2 : 0];

	// Copy the stats and geometries to the new shapefile for
	// each Stat.
	OGRPolygon *poly;
	for(auto it = stats.begin(); it != stats.end(); ++it) {
		feat = OGRFeature::CreateFeature(layer->GetLayerDefn());
		// Set original fields first.
		for(int i = 0; i < it->feat()->GetFieldCount(); ++i)
			feat->SetField(i, it->feat()->GetRawFieldRef(i));
		feat->SetField("stddev", it->stddev());
		feat->SetField("median", it->median());
		feat->SetField("sum", it->sum());
		feat->SetField("mean", it->mean());
		feat->SetField("count", it->count());
		if(numQuantiles > 0) {
			it->quantiles(quantiles, numQuantiles);
			for(int qt = 0; qt < numQuantiles + 2; ++qt)
				feat->SetField(("qt" + std::to_string(qt)).c_str(), quantiles[qt]);
		}
		poly = (OGRPolygon *) OGRGeometryFactory::createFromGEOS(gctx, (GEOSGeom) it->geom());
		feat->SetGeometry(poly);
		OGRErr err = layer->CreateFeature(feat);
		if(0 != err)
			throw "Failed to create feature.";
		layer->SetFeature(feat);
		OGRFeature::DestroyFeature(feat);
	}

}

void doRaster(std::string &raster, std::string &outfile, int band,
		std::vector<std::string> &files, std::vector<int> &classes, int numQuantiles, int outType) {

	GDALDataset *ds = (GDALDataset *) GDALOpen(raster.c_str(), GA_ReadOnly);
	if(ds == NULL)
		throw "Failed to open raster.";

	GDALRasterBand *b = ds->GetRasterBand(band);
	if(b == NULL)
		throw "Failed to load band.";

	int cols = ds->GetRasterXSize();
	int rows = ds->GetRasterYSize();
	double trans[6];
	ds->GetGeoTransform(trans);
	const char *proj = ds->GetProjectionRef();

	double minx = trans[0];
	double maxy = trans[3];
	double maxx = minx + (cols * trans[1]);
	double miny = maxy + (rows * trans[5]);

	std::map<int, std::map<int, Stat * > > stats;

	Grid<char> dg(cols, rows);

	if(0 != b->RasterIO(GF_Read, 0, 0, cols, rows, dg.grid(), cols, rows, GDT_Byte, 0, 0))
		throw "Failed to read from classification raster.";

	// Iterate over LAS files, process each one.
	las::ReaderFactory rf;
	std::cerr << "Processing files." << std::endl;
	for(unsigned int i = 0; i < files.size(); ++i) {

		std::ifstream in(files[i].c_str(), std::ios::in | std::ios::binary);
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();

		std::cout << "Processing file " << files[i]<< std::endl;

		int cls;
		while(r.ReadNextPoint()) {
			las::Point pt = r.GetPoint();

			if(!Util::inList(classes, (cls = pt.GetClassification().GetClass())))
				continue;

			double x = pt.GetX();
			double y = pt.GetY();
			double z = pt.GetZ();

			if(x < minx || x > maxx || y < miny || y > maxy)
				continue;

			int c = (int) ((x - minx) / (maxx - minx) * cols);
			int r = (int) ((maxy - y) / (maxy - miny) * rows);
			int val = (int) dg(c, r);

			Stat *st = stats[val][cls];
			if(st == NULL)
				st = stats[val][cls] = new Stat(val, cls);
			st->add(z);
		}

		in.close();
	}

	if(outType == CSV) {

		std::stringstream head;
		head << "id,cls,count,sum,min,max,mean,median,variance,stddev";
		for(int i = 0; i < numQuantiles + 2; ++i)
			head << ",q" << i;
		head << std::endl;

		std::ofstream out;
		out.open(outfile.c_str());
		out << head.str();
		out << std::setprecision(9);

		double quantiles[numQuantiles > 0 ? numQuantiles + 2 : 0];

		for(std::map<int, std::map<int, Stat*> >::iterator it = stats.begin(); it != stats.end(); ++it) {
			for(unsigned int c = 0; c < classes.size(); ++c) {
				int id = it->first;
				int cls = classes[c];
				Stat *st = stats[id][cls];
				out << id << "," << cls << "," << st->count() << "," << st->sum() << ","
					<< st->min() << "," << st->max() << "," << st->mean() << ","
					<< st->median() << "," << st->variance() << "," << st->stddev();
				if(numQuantiles > 0) {
					st->quantiles(quantiles, numQuantiles);
					for(int i = 0; i < numQuantiles + 2; ++i)
						 out << "," << quantiles[i];
				}
				out << std::endl;
				delete st;
			}
		}

		out.close();

	} else if(outType == RASTER) {

		int numBands = Stat::numBands(numQuantiles);
		std::vector<std::string> bandNames(numBands);
		Stat::bandNames(bandNames, numQuantiles);
		std::vector<double> values(numBands);

		int numRasterBands = numBands * classes.size();
		GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset *dst = drv->Create(outfile.c_str(), cols, rows, numRasterBands, GDT_Float32, NULL);
		GDALRasterBand *band = NULL;

		dst->SetGeoTransform(trans);
		dst->SetProjection(proj);

		Grid<float> fg(cols, rows);

		for(auto it = stats.begin(); it != stats.end(); ++it) {
			for(unsigned int c = 0; c < classes.size(); ++c) {
				int id = it->first;
				int cls = classes[c];
				Stat *st = stats[id][cls];
				st->bandValues(values, numQuantiles);
				for(int b = 0; b < numBands; ++b) {
					band = dst->GetRasterBand((c * numBands) + (b + 1));
					if(band == NULL)
						throw "Failed to retrieve raster band.";
					if(0 != band->RasterIO(GF_Read, 0, 0, cols, rows, fg.grid(), cols, rows, GDT_Float32, 0, 0))
						throw "Failed to read band raster.";
					for(int i = 0; i < rows * cols; ++i) {
						if(id == dg[i])
							fg[i] = values[b];
					}
					if(0 != band->RasterIO(GF_Write, 0, 0, cols, rows, fg.grid(), cols, rows, GDT_Float32, 0, 0))
						throw "Failed to write band raster.";
				}
				delete st;
			}
		}
	}

}

// TODO: If shapefile not given, computes stats for the entire point cloud.
void usage() {
	std::cerr << "Produces a zonal stats style output for points within polygons or a raster." << std::endl;
	std::cerr << "Usage: lasstats [options] lasfiles*" << std::endl;
	std::cerr << " -o Speficy an output file. This can be a shapefile (shp), a raster (tif) or a csv (csv)." << std::endl;
	std::cerr << " -s Specify a shapefile. Stats are computed for each polygon." << std::endl;
	std::cerr << " -r Specify a raster file. Stats are computed for each unique cell ID." << std::endl;
	std::cerr << " -l Specify a layer (shapefile only.)" << std::endl;
	std::cerr << " -b Specify a band (raster only.)" << std::endl;
	std::cerr << " -c A comma-delimited list of classes. If not given, all " << std::endl;
	std::cerr << "    classes are computed. Stats for all classes are computed" << std::endl;
	std::cerr << "    separately in either case." << std::endl;
	std::cerr << " -q If specified, is the number of quantiles to compute. Default 4." << std::endl;
}

int main(int argc, char ** argv) {

	std::vector<std::string> files;
	std::string outfile;
	std::string infile;
	int inType = 0;
	int outType = 0;
	std::string layername;
	int band = 1;
	int numQuantiles = 4;
	std::set<int> classSet;

	for(int i=1;i<argc;++i) {
		std::string arg(argv[i]);
		if(arg == "-s") {
			infile.assign(argv[++i]);
			inType = VECTOR;
		} else if(arg == "-o") {
			outfile.assign(argv[++i]);
		} else if(arg == "-c") {
			Util::intSplit(classSet, argv[++i]);
		} else if(arg == "-l") {
			layername.assign(argv[++i]);
		} else if(arg == "-r") {
			infile.assign(argv[++i]);
			if(inType != 0)
				throw "Raster and vector are mutually exclusive.";
			inType = RASTER;
		} else if(arg == "-b") {
			band = atoi(argv[++i]);
		} else if(arg == "-q") {
			numQuantiles = atoi(argv[++i]);
		} else {
			files.push_back(argv[i]);
		}
	}

	std::vector<int> classes(classSet.size());
	std::copy(classSet.begin(), classSet.end(), classes.begin());
	std::sort(classes.begin(), classes.end());

	std::string tmpout(outfile);
	std::transform(tmpout.begin(), tmpout.end(), tmpout.begin(), ::tolower);
	if(alg::ends_with(tmpout, ".shp")) {
		outType = VECTOR;
	} else if(alg::ends_with(tmpout, ".csv")) {
		outType = CSV;
	} else if(alg::ends_with(tmpout, ".tif")) {
		outType = RASTER;
	} else {
		std::cerr << "There are three allowed output types: .shp, .tif and .csv." << std::endl;
		usage();
		return 1;
	}

	if(outfile.empty()) {
		std::cerr << "An output file (-o) is required." << std::endl;
		usage();
		return 1;
	}

	if((inType != RASTER && inType != VECTOR) || infile.empty()) {
		std::cerr << "A shape file (-s) or raster (-r) is required, but not both." << std::endl;
		usage();
		return 1;
	}

	if(files.size() == 0) {
		std::cerr << "At least one input file is required." << std::endl;
		usage();
		return 1;
	}

	if(classes.size() == 0) {
		std::cerr << "WARNING! No classes specified. Matching all classes!" << std::endl;
	}

	if(band < 1) {
		std::cerr << "A band >= 1 is required." << std::endl;
		usage();
		return 1;
	}

	if(numQuantiles < 2) {
		std::cerr << "The number of quantiles must be at least 2." << std::endl;
		usage();
		return 1;
	}

	// We want the number of *dividers*.
	numQuantiles--;

	GDALAllRegister();

	int ret = 1;
	try {
		if(inType == RASTER) {
			doRaster(infile, outfile, band, files, classes, numQuantiles, outType);
		} else {
			doVector(infile, outfile, layername, files, classes, numQuantiles, outType);
		}
		ret = 0;
	} catch(const char *e) {
		std::cerr << e << std::endl;
		usage();
	}

	return ret;
}
