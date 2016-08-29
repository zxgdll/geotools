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

class Px {
public:
	double x;
	double y;
	std::vector<unsigned int> dn;
	
	Px(double x, double y) :
		x(x), y(y) {}

	void print(std::ostream &out, unsigned int id) {
		out << id << "," << x << "," << y;
		for(const unsigned int &v : dn)
			out << "," << v;
		out << std::endl;
	}

	bool isNoData(unsigned int nodata) {
		return dn[0] == nodata;
	}
};

class Poly {
public:
	unsigned int id;
	unsigned int nodata;
	std::unordered_map<size_t, std::unique_ptr<Px> > px;

	Poly(unsigned int id, unsigned int nodata) :
		id(id),
		nodata(nodata) {}

	void add(int col, int row, double x, double y, 	unsigned int dn) {
		size_t idx = ((size_t) col << 32 | row);
		if(px.find(idx) == px.end())
			px[idx] = std::unique_ptr<Px>(new Px(x, y));
		px[idx]->dn.push_back(dn);
	}

	void print(std::ostream &out) {
		for(const auto &it : px) {
			if(!it.second->isNoData(nodata))
				it.second->print(out, id);
		}
	}
};

	
void processSpectralFile(const SpectralConfig &config, Raster<unsigned int> &idxRaster, const std::string &specFilename) {
	g_loglevel(G_LOG_DEBUG);

	std::set<int> bands;
	bands.insert(config.bands.begin(), config.bands.end());
	if(bands.empty()) {
		Raster<unsigned int> tmp(specFilename);
		for(int i = 1; i <= tmp.bandCount(); ++i)
			bands.emplace(i);
	}
	g_debug("bands " << bands.size());

	std::unordered_map<unsigned int, std::unique_ptr<Poly> > polys;
	std::unordered_set<unsigned int> keep;
	std::unordered_set<unsigned int> skip;

	for(int row = 0; row < idxRaster.rows(); ++row) {
		keep.clear();

		// Iterate over bands to populate Polys
		for(const int &band : bands) {

			Raster<unsigned short> specRaster(specFilename, band);

			for(int col = 0; col < idxRaster.cols(); ++col) {
				unsigned int id = idxRaster.get(col, row);
				if(id == 0)
					continue;
				keep.emplace(id);
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

