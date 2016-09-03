#include <unordered_set>
#include <unordered_map>
#include <list>
#include <memory>
#include <iostream>
#include <iomanip>

#include "geotools.h"
#include "raster.hpp"
#include "spectral.hpp"
#include "sqlite.hpp"

using namespace geotools::spectral;
using namespace geotools::spectral::config;
using namespace geotools::db;

#include "raster.cpp"

/**
 * Represents a single pixel across all bands.
 * Contains the x/y position of the pixel center, and a list
 * of digital numbers corresponding to the bands, in order.
 */
class Px {
public:
	double x;
	double y;
	unsigned int id;
	std::map<int, unsigned short> dn;
	
	/**
	 * Construct the Px with the position.
	 */
	Px(unsigned int id, double x, double y) :
		id(id), x(x), y(y) {}

	/**
	 * Print to the output stream. Prints the x, y and
	 * all band in order, comma-delimited.
	 */
	void print(std::ostream &out) {
		out << id << "," << x << "," << y;
		for(const auto &it : dn)
			out << "," << it.second;
		out << std::endl;
	}

	void save(SQLite &db) {
		std::map<std::string, std::string> fields;
		fields["id"] = std::to_string(id);
		for(const auto &it : dn)
			fields[std::to_string(it.first)] = std::to_string(it.second);
		db.addPoint(x, y, 0, fields);
	}

	/**
	 * Returns true if the pixek should be abandoned because it contains nodata.
	 */
	bool isNoData(unsigned int nodata) {
		// TODO: Consult on the best way to determine whether a pixel is
		// abandoned because of nodata.
		return dn[0] == nodata;
	}
};

/**
 * Computes a band list. If the given list is empty, populates the output list with
 * all of the avaible bands. Otherwise returns a list of the desired bands.
 */
std::vector<int> computeBandList(const std::set<int> &bands, const std::string &specFilename) {
	std::vector<int> b;
	b.assign(bands.begin(), bands.end());
	if(b.empty()) {
		Raster<unsigned int> tmp(specFilename);
		for(int i = 1; i <= tmp.bandCount(); ++i)
			b.push_back(i);
		g_warn("No bands selected; using all " << b.size() << " bands.");
	}
	return b;
}

/**
 * Computes the bounding box which describes the overlapping area between the index
 * raster and the spectral raster.
 */
Bounds computeOverlapBounds(Raster<unsigned int> &idxRaster, const std::string &specFilename) {
	Raster<float> tmp(specFilename);
	return idxRaster.bounds().intersection(tmp.bounds());
}

/**
 * Process a single spectral file.
 * config       -- A SpectralConfig object with operational parameters.
 * idxRaster    -- A Raster containing the polygon IDs that are used to group the extracted values.
 * specFilename -- The filename of the current spectral file.
 */
void processSpectralFile(const SpectralConfig &config, const std::vector<int> &bands, const std::string &specFilename, SQLite &db) {

	g_debug("processSpectralFile [config] [raster] " << specFilename);

	int startRow, endRow, startCol, endCol;
	unsigned int idxNodata = config.idxNodata;
	unsigned short specNodata = config.specNodata;
	Bounds bounds;

	{
		Raster<unsigned int> idxRaster(config.indexFilename);
		// Get the start and end cols/rows from the bounds intersection.
		bounds.collapse();
		bounds.extend(computeOverlapBounds(idxRaster, specFilename));
		startRow = idxRaster.toRow(idxRaster.resolutionY() > 0 ? bounds.miny() : bounds.maxy());
		endRow = idxRaster.toRow(idxRaster.resolutionY() > 0 ? bounds.maxy() : bounds.miny());
		startCol = idxRaster.toCol(bounds.minx());
		endCol = idxRaster.toCol(bounds.maxx());
		// If nodata is not set, get it from the raster.
		if(!config.hasIdxNodata)
			idxNodata = idxRaster.nodata();
	}

	// If nodata is not set, get it from the raster.
	if(!config.hasSpecNodata) {
		Raster<unsigned short> specRaster(specFilename);
		specNodata = specRaster.nodata();
	}

	g_debug(" - rows: " << startRow << " -> " << endRow << "; cols: " << startCol << " -> " << endCol);
	g_debug(" - bands " << bands.size() << "; bounds: " << bounds.print());

	//std::unordered_map<unsigned int, std::unique_ptr<Poly> > polys;
	std::unordered_map<size_t, std::unique_ptr<Px> > px;
	std::unordered_set<size_t> skip;
	std::unordered_set<size_t> finalize;

	// Iterate over bands to populate Polys
	#pragma omp parallel for
	for(int block = startRow; block < endRow; block += 100) {

		Raster<unsigned int> idxRaster(config.indexFilename);
		Raster<unsigned short> specRaster(specFilename, bands[0]);
		
		for(int row = block; row < g_min(endRow, block + 100); ++row) {
			for(int b = 0; b < bands.size(); ++b) {
				int band = bands[b];
				specRaster.setBand(band);
				for(int col = startCol; col < endCol; ++col) {
					unsigned int id = idxRaster.get(col, row);
					if(id == idxNodata)
						continue;
					double x = idxRaster.toX(col) + idxRaster.resolutionX() / 2.0;
					double y = idxRaster.toY(row) + idxRaster.resolutionY() / 2.0;
					int scol = specRaster.toCol(x);
					int srow = specRaster.toRow(y);
					if(specRaster.has(scol, srow)) {
						unsigned int v = specRaster.get(scol, srow);
						if(v == specNodata)
							continue;
						size_t idx = ((size_t) col << 32) | row;
						#pragma omp critical (PX)
						{
							if(skip.find(idx) == skip.end()) {
								if(px.find(idx) == px.end())
									px[idx] = std::unique_ptr<Px>(new Px(id, x, y));
								Px &p = *(px[idx]);
								p.dn[band] = v;
								//g_debug("px id: " << id << "; idx: " << idx << "; dn size: " << p.dn.size());
								if(p.dn.size() == bands.size())
									finalize.insert(idx);
							} else {
								skip.emplace(idx);
							}
						}
					}
				}
			}
		}
		// Remove polys that have IDs that are not in the current row.	
		std::cout << std::fixed << std::setprecision(3);
		#pragma omp critical (PX)
		{
			g_debug("Save " << finalize.size());
			db.begin();
			for(const size_t &idx : finalize) {
				if(px[idx]->dn.size() == bands.size()) 
					px[idx]->save(db);
			}
			db.commit();
		}
		#pragma omp critical (PX)
		{
			g_debug("Finalize.");
			for(const size_t &idx : finalize)
				px.erase(idx);
			finalize.clear();
		}

	}
}

void Spectral::extractSpectra(const SpectralConfig &config) {
	g_debug("extractSpectra [config]");

	// Check the config for problems.
	config.check();
	
	// Check the bands list; fill it if necessary.
	std::vector<int> bands = computeBandList(config.bands, config.spectralFilenames[0]);
	
	// Develop the fields list for the DB.
	std::map<std::string, int> fields;
	fields["id"] = SQLite::INTEGER;
	for(const int &band : bands)
		fields["b" + std::to_string(band)] = SQLite::INTEGER;

	// Initialize the SQLite table.
	SQLite db(config.outputFilename, SQLite::POINT, config.srid, fields);

	// Do the work.
	for(const std::string &specFilename : config.spectralFilenames) {
		g_debug(" - file: " << specFilename);
		processSpectralFile(config, bands, specFilename, db);
	}
}

