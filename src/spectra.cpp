#include <unordered_set>
#include <unordered_map>
#include <list>
#include <memory>
#include <iostream>
#include <iomanip>

#include "geotools.h"
#include "raster.hpp"

#include "spectra.hpp"

using namespace geotools::spectral;
using namespace geotools::spectral::config;

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
	std::vector<unsigned int> dn;
	
	/**
	 * Construct the Px with the position.
	 */
	Px(double x, double y) :
		x(x), y(y) {}

	/**
	 * Print to the output stream. Prints the x, y and
	 * all band in order, comma-delimited.
	 */
	void print(std::ostream &out) {
		out << x << "," << y;
		for(const unsigned int &v : dn)
			out << "," << v;
		out << std::endl;
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
 * Represents a contiguous polygon comprised of pixels having the same ID.
 * Contains Px object (indexed by a pixel index) which in turn contain
 * spectral band values.
 */
class Poly {
public:
	unsigned int id;
	unsigned int nodata;
	std::unordered_map<size_t, std::unique_ptr<Px> > px;

	/**
	 * Construct a poly with the given polygon ID and nodata value.
	 */
	Poly(unsigned int id, unsigned int nodata) :
		id(id),
		nodata(nodata) {}

	/**
	 * Add a band value at the given column and row. The x/y position is also given.
	 */
	void add(int col, int row, double x, double y, 	unsigned int dn) {
		size_t idx = ((size_t) col << 32 | row);
		if(px.find(idx) == px.end())
			px[idx] = std::unique_ptr<Px>(new Px(x, y));
		px[idx]->dn.push_back(dn);
	}

	/**
	 * Print out the polygons spectral values for each pixel in CSV format.
	 */
	void print(std::ostream &out) {
		for(const auto &it : px) {
			if(!it.second->isNoData(nodata)) {
				out << id << ",";
				it.second->print(out);
			}
		}
	}
};

std::set<int> computeBandList(const std::set<int> &bands, const std::string &specFilename) {
	std::set<int> b;
	b.insert(bands.begin(), bands.end());
	if(b.empty()) {
		Raster<unsigned int> tmp(specFilename);
		for(int i = 1; i <= tmp.bandCount(); ++i)
			b.emplace(i);
	}
	return b;
}

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
void processSpectralFile(const SpectralConfig &config, Raster<unsigned int> &idxRaster, const std::string &specFilename) {
	g_loglevel(G_LOG_DEBUG);
	g_debug("processSpectralFile [config] [raster] " << specFilename);

	// TODO: More versatile band selection.
	std::set<int> bands = computeBandList(config.bands, specFilename);
	Bounds bounds = computeOverlapBounds(idxRaster, specFilename);

	g_debug(" - bands " << bands.size() << "; bounds: " << bounds.print());

	std::unordered_map<unsigned int, std::unique_ptr<Poly> > polys;
	std::unordered_set<unsigned int> keep;
	std::unordered_set<unsigned int> skip;

	int startRow = idxRaster.toRow(bounds.miny());
	int endRow = idxRaster.toRow(bounds.maxy());
	int startCol = idxRaster.toCol(bounds.minx());
	int endCol = idxRaster.toCol(bounds.maxx());

	g_debug(" - rows: " << startRow << " -> " << endRow << "; cols: " << startCol << " -> " << endCol);

	for(int row = startRow; row < endRow; ++row) {
		keep.clear();

		// Iterate over bands to populate Polys
		for(const int &band : bands) {

			Raster<unsigned short> specRaster(specFilename, band);

			for(int col = startCol; col < endCol; ++col) {
				unsigned int id = idxRaster.get(col, row);
				if(id == idxRaster.nodata())
					continue;
				keep.emplace(id);						// Keeps track of IDs found in the row; if a poly isn't in this list, output it.
				double x = idxRaster.toX(col) + idxRaster.resolutionX() / 2.0;
				double y = idxRaster.toY(row) + idxRaster.resolutionY() / 2.0;
				int scol = specRaster.toCol(x);
				int srow = specRaster.toRow(y);
				if(specRaster.has(scol, srow)) {
					unsigned int v = specRaster.get(scol, srow);
					if(v == config.nodata) {
						polys.erase(id);
						skip.emplace(id);
					} else if(skip.find(id) == skip.end()) {
						if(polys.find(id) == polys.end())
							polys[id] = std::unique_ptr<Poly>(new Poly(id, config.nodata));
						polys[id]->add(col, row, x, y, specRaster.get(scol, srow));
					}
				}
			}
		}
		
		// Remove polys that have IDs that are not in the current row.	
		std::list<unsigned int> premove;
		for(const auto &it : polys) {
			if(keep.find(it.first) == keep.end()) {
				it.second->print(std::cout << std::fixed << std::setprecision(12));
				premove.push_back(it.first);
			}
		}
		for(const unsigned int &it : premove)
			polys.erase(it);
	}
}

void Spectral::extractSpectra(const SpectralConfig &config) {
	Raster<unsigned int> idxRaster(config.indexFilename);
	for(const std::string &specFilename : config.spectralFilenames) {
		processSpectralFile(config, idxRaster, specFilename);
	}
}

