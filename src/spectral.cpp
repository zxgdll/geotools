#include <unordered_set>
#include <unordered_map>
#include <list>
#include <memory>
#include <iostream>
#include <iomanip>

#include "geotools.h"
#include "raster.hpp"

#include "spectral.hpp"

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

	/**
	 * Return the number of band DNs in the px object.
	 */
	int count() {
		const Px *p = px.begin()->second.get();
		return p->dn.size();
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
void processSpectralFile(const SpectralConfig &config, const std::string &specFilename) {

	g_debug("processSpectralFile [config] [raster] " << specFilename);

	// TODO: More versatile band selection.
	std::vector<int> bands = computeBandList(config.bands, specFilename);
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

	// Iterate over bands to populate Polys
	#pragma omp parallel for
	for(int b = 0; b < bands.size(); ++b) {

		int band = bands[b];

		Raster<unsigned int> idxRaster(config.indexFilename);
		Raster<unsigned short> specRaster(specFilename, band);
		std::unordered_set<unsigned int> keep;
		std::unordered_set<unsigned int> skip;
		std::unordered_map<unsigned int, std::unique_ptr<Poly> > polys;

		for(int row = startRow; row < endRow; ++row) {
			keep.clear();

			for(int col = startCol; col < endCol; ++col) {
				unsigned int id = 0;
				#pragma omp critical(idRead)
				{ 
					id = idxRaster.get(col, row); // TODO: Use a raster for each thread.
				}
				if(id == idxNodata)
					continue;
				keep.emplace(id); // Keeps track of IDs found in the row; if a poly isn't in this list, output it.
				double x = idxRaster.toX(col) + idxRaster.resolutionX() / 2.0;
				double y = idxRaster.toY(row) + idxRaster.resolutionY() / 2.0;
				int scol = specRaster.toCol(x);
				int srow = specRaster.toRow(y);
				if(specRaster.has(scol, srow)) {
					unsigned int v = specRaster.get(scol, srow);
					if(v == specNodata) {
						polys.erase(id);
						skip.emplace(id);
					} else if(skip.find(id) == skip.end()) {
						if(polys.find(id) == polys.end())
							polys[id] = std::unique_ptr<Poly>(new Poly(id, specNodata));
						polys[id]->add(col, row, x, y, specRaster.get(scol, srow));
					}
				}
			}
			// Remove polys that have IDs that are not in the current row.	
			std::list<unsigned int> premove;
			for(const auto &it : polys) {
				if(keep.find(it.first) == keep.end() && it.second->count() == bands.size()) {
					#pragma omp critical(output)
					{
						it.second->print(std::cout << std::fixed << std::setprecision(3));
					}
					premove.push_back(it.first);
				}
			}
			for(const unsigned int &it : premove)
				polys.erase(it);
		}
	}
}

void Spectral::extractSpectra(const SpectralConfig &config) {
	g_debug("extractSpectra [config]");

	config.check();

	for(const std::string &specFilename : config.spectralFilenames) {
		g_debug(" - file: " << specFilename);
		processSpectralFile(config, specFilename);
	}
}

