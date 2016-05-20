#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <liblas/liblas.hpp>
#include <proj_api.h>
#include <gdal_priv.h>
#include <gdal.h>
#include <cpl_conv.h> // for CPLMalloc()

#include "geotools.h"
#include "Util.hpp"
#include "Raster.hpp"

#define LAS2CSRS_DATA "LAS2CSRS_DATA"

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;
namespace las = liblas;

bool quiet = true;

// Binary interpolation.
double _binterp(float *grid, double c, double r, int c0, int r0, int c1, int r1, int width) {
	double x1 = (c1 - c) / (c1 - c0) * grid[r0 * width + c0] + (c - c0) / (c1 - c0) * grid[r0 * width + c1];
	double x2 = (c1 - c) / (c1 - c0) * grid[r1 * width + c0] + (c - c0) / (c1 - c0) * grid[r1 * width + c1];
	return (r1 - r) / (r1 - r0) * x1 + (r - r0) / (r1 - r0) * x2;
}

/**
 * Convert projected distance in mm to latlon in radians.
 * dx, dy, dx  -- Distance in m.
 * lat         -- The latitude at which distances are computed
 * h           -- Ellipsoidal height.
 * a           -- Semi-major axis
 * e2          -- Eccentricity^2
 * count       -- Number of points.
 * dlat        -- The distance in rad (out).
 * dlon        -- The distance in rad (out).
 */
void _shift2latlon(Grid<double> &gdx, Grid<double> &gdy, Grid<double> &glat, Grid<double> &gh,
	double a, double e2, int count, Grid<double> &gdlat, Grid<double> &gdlon) {

	double *dx = gdx.grid();
	double *dy = gdy.grid();
	double *lat = glat.grid();
	double *dlat = gdlat.grid();
	double *dlon = gdlon.grid();
	double *h = gh.grid();

	double r, l, m, n;
	for(int i = 0; i < count; ++i) {
		l = *(lat + i);
		m = a * (1.0 - e2) / pow((1.0 - e2 * _sq(sin(l))), 3.0/2.0); 	// Meridional radius of curvature.
		n = a / pow((1.0 - e2 * _sq(sin(l))), 1.0/2.0); 				// Parallel radius of curvature.
		r = n * cos(l); 												// Radius of parallel.
		*(dlon + i) = *(dx + i) / (r + *(h + i));
		*(dlat + i) = *(dy + i) / (m + *(h + i));
	}
}

/**
 * Loads and interpolates the NAD83(CSRS) shift grid.
 */
class ShiftGrid {
public:

	/**
	 * Load the shift grid raster and store the x, y, z shifts internally.
	 */
	 // TODO: Only load the portion of the grid necessary for the point cloud.
	 //       On the other hand, the grid is so small, who cares?
	void load() {

		char path[PATH_MAX];
		sprintf(path, "%s%s", std::getenv("LAS2CSRS_DATA"), "/NAD83v6VG.tif");

		GDALAllRegister();
		GDALDataset *ds;
		GDALRasterBand *xband, *yband, *zband;

		ds = (GDALDataset *) GDALOpen(path, GA_ReadOnly);
		if(!ds)
			throw "Failed to load shift grid.";

		xband = ds->GetRasterBand(1);
		yband = ds->GetRasterBand(2);
		zband = ds->GetRasterBand(3);
		if(!xband || !yband || !zband)
			throw "Failed to retrieve shift bands.";

		width = xband->GetXSize();
		height = yband->GetYSize();
		if(width <= 0 || height <= 0)
			throw "The dimensions of the shift grid are invalid.";

		xg.init(width, height);
		yg.init(width, height);
		zg.init(width, height);

		CPLErr xe = xband->RasterIO(GF_Read, 0, 0, width, height, xg.grid(), width, height, GDT_Float32, 0, 0);
		CPLErr ye = yband->RasterIO(GF_Read, 0, 0, width, height, yg.grid(), width, height, GDT_Float32, 0, 0);
		CPLErr ze = zband->RasterIO(GF_Read, 0, 0, width, height, zg.grid(), width, height, GDT_Float32, 0, 0);

		if(xe == CE_Failure || ye == CE_Failure || ze == CE_Failure)
			throw "Failed to read shift grid.";

		tg.init(6, 1);
		ds->GetGeoTransform(tg.grid());
	}

	/**
	 * Compute the the shifts in x, y and z at position x, y. x and y are radians,
	 * dx, dy and dz are m.
	 */
	void interpolate(Grid<double> &gx, Grid<double> &gy, Grid<double> &gdx, 
		Grid<double> &gdy, Grid<double> &gdz, int count) {

		double *x = gx.grid();
		double *y = gy.grid();
		double *dx = gdx.grid();
		double *dy = gdy.grid();
		double *dz = gdz.grid();

		double c, r;
		int c0, c1, r0, r1;
		for(;count;--count, ++x, ++y, ++dx, ++dy, ++dz) {
			c = ((_deg(*x) - tg[0]) / tg[1]);
			r = ((_deg(*y) - tg[3]) / tg[5]);
			c0 = (int) c;
			r0 = (int) r;
			c1 = c0 + 1;
			r1 = r0 + 1;
			if(c0 < 0) c0 = 0;
			if(r0 < 0) r0 = 0;
			if(c1 >= width) c1 = width-1;
			if(r1 >= height) r1 = height-1;
			*dx = _binterp(xg.grid(), c, r, c0, r0, c1, r1, width) / 1000.0;
			*dy = _binterp(yg.grid(), c, r, c0, r0, c1, r1, width) / 1000.0;
			*dz = _binterp(zg.grid(), c, r, c0, r0, c1, r1, width) / 1000.0;
		}
	}

	~ShiftGrid() {
	}

private:
	MemRaster<float> xg;
	MemRaster<float> yg;
	MemRaster<float> zg;
	MemRaster<double> tg;
	int width;
	int height;

};

class Params {
public:
	// Source reference frame.
	std::string ffrom;
	// From and to epochs.
	double efrom;
	double eto;
	// SRIDs of the from and to CRSes.
	int fsrid;
	int tsrid;
	// Transform parameters; loaded from the itrf file.
	double tx, ty, tz, dtx, dty, dtz;		// Shifts, rates.
	double rx, ry, rz, drx, dry, drz;		// Rotations, rates.
	double epoch;							// ITRF Transform epoch.
	//double dt; 								// Time delta
	double d, dd; 							// Scale, scale rate.

};

/**
 * Performs the work of transforming coordinates from a given reference frame to NAD83(CSRS).
 */
class Transformer {
public:

	/**
	 * Prepares the Helmert transformation parameters for the given transformation.
	 * These will be used in later method calls.
	 * The constructor will load an process the grid shift file and the transformation database.
	 * ffrom -- The name of the reference frame, e.g. 'itrf90'
	 * efrom -- The epoch of data collection (decimal years), e.g. 1994.2
	 * eto   -- The target epoch. One might select 1997.0 for BC or 2002.0 for Albera (etc.)
	 * fsrid -- The SRID (EPSG code) of the source.
	 * tsrid -- The SRID of the destination. This is the code for the UTM zone under
	 *          NAD83(CSRS), e.g. 2956 for UTM12N.
	 */
	Transformer(std::string &ffrom, float efrom, float eto, int fsrid, int tsrid) {

		if(!quiet) {
			std::cerr << "Transformer:" << std::endl
				<< " -- src ref frame: " << ffrom << "; src epoch: " << efrom << "; dst epoch: " << eto << std::endl
				<< " -- src srid: " << fsrid << "; dst srid: " << tsrid << std::endl;
		}
		params.efrom = efrom;
		params.eto = eto;
		params.fsrid = fsrid;
		params.tsrid = tsrid;
		params.ffrom.assign(ffrom);

		initProjections();
		loadHelmert(params);
		shiftGrid.load();
	}

	/**
	 * Transforms coordinate(s) from one available reference frame to NAD83(CSRS).
	 * x, y, z -- Coordinate arrays.
	 * count   -- The number of coordinates.
	 * bounds  -- The bounds of the transformed coordinate list.
	 */
	void transformPoints(Grid<double> &x, Grid<double> &y, Grid<double> &z, int count, double bounds[6]) {

		if(!quiet) {
			std::cerr << std::setprecision(9) << "transformPoints" << std::endl;
			std::cerr << "  -- Original: " << x.grid()[0] << ", " << y.grid()[0] << ", " << z.grid()[0] << std::endl;
		}

		// Copy the coordinate arrays for transformation.
		MemRaster<double> x0(count, 1);
		MemRaster<double> y0(count, 1);
		MemRaster<double> z0(count, 1);
		memcpy(x0.grid(), x.grid(), sizeof(double) * count);
		memcpy(y0.grid(), y.grid(), sizeof(double) * count);
		memcpy(z0.grid(), z.grid(), sizeof(double) * count);

		// Project to Cartesian 3D. (c)
		pj_transform(p1, p2, count, 1, x.grid(), y.grid(), z.grid());

		// Transform to latlon. (b)
		pj_transform(p1, p4, count, 1, x0.grid(), y0.grid(), z0.grid());

		if(!quiet)
			std::cerr << "  -- Cartesian: " << x.grid()[0] << ", " << y.grid()[0] << ", " << z.grid()[0] << std::endl;

		if(!quiet)
			std::cerr << "  -- Latlon: " << _deg(x0.grid()[0]) << ", " << _deg(y0.grid()[0]) << ", " << z0.grid()[0] << std::endl;

		// Transform to csrs using Helmert (etc.) params. (c)
		epochTransform(x, y, z, count, params.efrom - 1997.0);

		if(!quiet)
			std::cerr << "  -- Cartesian + epoch: " << x.grid()[0] << ", " << y.grid()[0] << ", " << z.grid()[0] << std::endl;

		// Only use the grid shift if the epoch changes.
		if(params.efrom != params.eto) {

			// Initalize shift arrays.
			MemRaster<double> dx(count, 1);
			MemRaster<double> dy(count, 1);
			MemRaster<double> dz(count, 1);
			
			// Interpolate shifts using latlon coords -- returns in mm. (d)
			shiftGrid.interpolate(x0, y0, dx, dy, dz, count);
			
			if(!quiet)
				std::cerr << "  -- Velocities: " << dx.grid()[0] << ", " << dy.grid()[0] << ", " << dz.grid()[0] << std::endl;

			// Transform mm shifts to latlon
			MemRaster<double> dlat(count, 1);
			MemRaster<double> dlon(count, 1);

			// Get projection's spheroid props.
			double a, e2;
			pj_get_spheroid_defn(p4, &a, &e2);

			// Apply shifts (mm) to angular coords.
			_shift2latlon(dx, dy, y0, z0, a, e2, count, dlat, dlon);

			if(!quiet)
				std::cerr << "  -- Lat lon shifts: " << dlat.grid()[0] << ", " << dlon.grid()[0] << std::endl;
			
			double dt = 1997.0 - params.eto;
			// Apply shifts to latlon coords.
			for(int i = 0; i < count; ++i) {
				*(x0.grid() + i) += *(dlon.grid() + i) * dt;
				*(y0.grid() + i) += *(dlat.grid() + i) * dt;
				*(z0.grid() + i) += *(dz.grid() + i) * dt;
			}
			
			if(!quiet)
				std::cerr << "  -- Lat lon shifted: " << _deg(x0.grid()[0]) << ", " << _deg(y0.grid()[0]) << ", " << _deg(z0.grid()[0]) << std::endl;

			// Transform latlon to target proj
			pj_transform(p4, p3, count, 1, x0.grid(), y0.grid(), z0.grid());
			
			if(!quiet)
				std::cerr << "  -- Final: " << x0.grid()[0] << ", " << y0.grid()[0] << ", " << z0.grid()[0] << std::endl;

			// Assign the shifted coords to the output arrays.
			memcpy(x.grid(), x0.grid(), sizeof(double) * count);
			memcpy(y.grid(), y0.grid(), sizeof(double) * count);
			memcpy(z.grid(), z0.grid(), sizeof(double) * count);

		} else {

			// Reproject to the dest coordinates
			pj_transform(p2, p3, count, 1, x.grid(), y.grid(), z.grid());

		}

		// Expand the bounds for the new header.
		for(int i = 0; i < count; ++i) {
			if(x[i] < bounds[0]) bounds[0] = x[i];
			if(x[i] > bounds[1]) bounds[1] = x[i];
			if(y[i] < bounds[2]) bounds[2] = y[i];
			if(y[i] > bounds[3]) bounds[3] = y[i];
			if(z[i] < bounds[4]) bounds[4] = z[i];
			if(z[i] > bounds[5]) bounds[5] = z[i];
		}			
	}

	/**
	 * Transform the given LAS file(s) from the given reference frame to the one
	 * configured in this Transformer. 
	 * srcfile  -- The folder containing las files, or the path to a single las file.
	 * dstdir   -- The destination folder.
	 */
	int transformLas(std::string &srcfile, std::string &dstdir, bool overwrite) {
		
		if(!quiet)
			std::cerr << "transformLas" << std::endl;

		const fs::path src(srcfile);
		const fs::path dst(dstdir);

		// Check the files for sanity.
		if(src == dst)
			throw "Destination and source are the same: " + dst.string() + "==" + src.string();

		if(!fs::exists(dst)) {
			if(!fs::create_directory(dst))
				throw "Failed to create output directory: " + dst.string();
		}

		// Get the list of las files.
		if(!quiet)
			std::cerr << " -- Getting file list." << std::endl;

		std::list<fs::path> files;
		if(fs::is_regular_file(src)) {
			if(overwrite || !fs::exists(dst / src.leaf()))
				files.push_back(src);
		} else {
			std::string ext(".las");
			fs::directory_iterator end;
			fs::directory_iterator di(src);
			for(; di != end; ++di) {
				std::string p(di->path().string());
				alg::to_lower(p);
				if((overwrite || !fs::exists(dst / di->path().leaf())) && alg::ends_with(p, ext)) {
					files.push_back(di->path());
				} else {
					throw "The destination file exists and -o is not specified.";
				}
			}
		}

		if(files.size() == 0)
			throw "No matching files found.";

		if(!quiet)
			std::cerr << "Processing " << files.size() << " files." << std::endl;

		// Start
		las::WriterFactory wf;
		las::ReaderFactory rf;

		// The overall bounds: min x, max x, min y, max y, min z, max z
		double bounds[] = { DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_POS, DBL_MAX_NEG };

		for(std::list<fs::path>::iterator it = files.begin(); it != files.end(); ++it) {

			const char * filename = it->c_str();

			std::cout << "Processing file " << filename << std::endl;

			// Open the source file.
			std::ifstream in(filename, std::ios::in | std::ios::binary);
			las::Reader r = rf.CreateWithStream(in);
			las::Header h = r.GetHeader();

			// Open the destination file.
			fs::path outfile(dst / it->leaf());
			las::Header dsth(h);

			int count = h.GetPointRecordsCount();

			MemRaster<double> x(count, 1);
			MemRaster<double> y(count, 1);
			MemRaster<double> z(count, 1);

			// Iterate over the points.
			for(int i = 0; r.ReadNextPoint(); ++i) {
				las::Point pt = r.GetPoint();
				*(x.grid() + i) = pt.GetX();
				*(y.grid() + i) = pt.GetY();
				*(z.grid() + i) = pt.GetZ();
			}

			// Transform the points in-place.
			transformPoints(x, y, z, count, bounds);

			// Set bounds.
			dsth.SetMin(bounds[0], bounds[2], bounds[4]);
			dsth.SetMax(bounds[1], bounds[3], bounds[5]);
			std::ofstream out(outfile.c_str(), std::ios::out | std::ios::binary);
			las::Writer w(out, dsth);

			// Iterate over the points.
			r.Reset();
			for(int i = 0; r.ReadNextPoint(); ++i) {
				las::Point pt = r.GetPoint();
				pt.SetX(*(x.grid() + i));
				pt.SetY(*(y.grid() + i));
				pt.SetZ(*(z.grid() + i));
				// Write to the output
				w.WritePoint(pt);
			}

			// Close files.
			in.close();
			out.close();

		}

		return 0;

	}

	~Transformer() {
		pj_free(p1);
		pj_free(p2);
		pj_free(p3);
		pj_free(p4);
	}

private:

	// Projection objects.
	projPJ p1;
	projPJ p2;
	projPJ p3;
	projPJ p4;

	// The shift grid.
	ShiftGrid shiftGrid;

	Params params;

	inline double _sec2rad(double x) {
		return (x * 0.000290888);
	}

	/**
		Transform the coordinate using the procedure listed in Craymer (2006).

		x, y, z -- 	The coordinate arrays.
		count 	-- 	The number of points.
	*/
	void epochTransform(Grid<double> &gx, Grid<double> &gy, Grid<double> &gz, int count, double dt) {

		if(!quiet)
			std::cerr << "epochTransform" << std::endl;

		double *x = gx.grid();
		double *y = gy.grid();
		double *z = gz.grid();

		double cx = params.tx + params.dtx * dt;            // Translation in X plus X velocity * time
		double cy = params.ty + params.dty * dt;
		double cz = params.tz + params.dtz * dt;
		double s = 1.0 + params.d + params.dd * dt;         // Scale plus delta scale * time
		double rx = _sec2rad(params.rx + params.drx * dt); // Rotation in X plus velocity * time; seconds to radians.
		double ry = _sec2rad(params.ry + params.dry * dt);
		double rz = _sec2rad(params.rz + params.drz * dt);

		if(!quiet) {
			std::cerr << " -- cx: " << cx << "; cy: " << cy << "; cz: " << cz << std::endl;
			std::cerr << " -- rx: " << rx << "; ry: " << ry << "; rz: " << rz << std::endl;
			std::cerr << " -- s: " << s << std::endl;
		}

		for(;count;--count, ++x, ++y, ++z) {
			*x = cx + s * (*x - rz * *y + ry * *z);
			*y = cy + s * (*y + rz * *x - rx * *z);
			*z = cz + s * (*z - ry * *x + rx * *y);
		}
	}

	/**
	 * Initialize the proj4 projection objects.
	 */
	void initProjections() {

		if(!quiet)
			std::cerr << "initProjection" << std::endl;

		// Initialize projections.
		if(!params.fsrid || !params.tsrid)
			throw "SRIDs are not set.";
		char str[128];
		sprintf(str, "+init=epsg:%u", params.fsrid);
		p1 = pj_init_plus(str);
		if(!p1)
			throw std::string(pj_strerrno(pj_errno));
		sprintf(str, "+init=epsg:%u", params.tsrid);
		p3 = pj_init_plus(str);
		if(!p3)
			throw std::string(pj_strerrno(pj_errno));
		p2 = pj_init_plus("+init=epsg:4978");
		if(!p2)
			throw std::string(pj_strerrno(pj_errno));
		p4 = pj_init_plus("+init=epsg:4326");
		if(!p4)
			throw std::string(pj_strerrno(pj_errno));
	}

	/**
	 * Load the transformation database.
	 */
	void loadHelmert(Params &p) {

		if(!quiet)
			std::cerr << "loadHelmert" << std::endl;

		char path[PATH_MAX];
		sprintf(path, "%s%s", std::getenv(LAS2CSRS_DATA), "/itrf.csv");

		char ffrom[64], fto[64];
		bool found = false;
		float epoch, tx, ty, tz, rx, ry, rz, d, dtx, dty, dtz, drx, dry, drz, dd;
		FILE *f = fopen(path, "r");
		if(f == NULL)
			throw "ITRF database file not found.";
		while(fscanf(f, " %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ",
						ffrom, fto, &epoch, &tx, &ty, &tz, &rx, &ry, &rz, &d,
						&dtx, &dty, &dtz, &drx, &dry, &drz, &dd) > 0) {
			std::string _ffrom(ffrom);
			if(!quiet)
				std::cerr << " -- Checking: " << _ffrom << std::endl;
			if(_ffrom == p.ffrom) {
				if(!quiet)
					std::cerr << " -- Found entry for " << _ffrom << std::endl;
				p.epoch = epoch;
				p.d = d / 1000000000.0; 		// Listed as ppb.
				p.dd = dd / 1000000000.0;
				p.tx = tx;
				p.ty = ty;
				p.tz = tz;
				p.rx = rx;
				p.ry = ry;
				p.rz = rz;
				p.dtx = dtx;
				p.dty = dty;
				p.dtz = dtz;
				p.drx = drx;
				p.dry = dry;
				p.drz = drz;
				found = true;
				break;
			}
		}
		fclose(f);

		if(!found)
			throw "Failed to find a transformation matching the parameters.";
	}
};

	
void usage() {
	std::cout << "Usage: las2csrs [options] <src file or dir> <dst dir> <src ref frame> <src epoch> <dst epoch> <srd srid> <dst srid>" << std::endl;
	std::cout << " -o     Overwrite existing files. Defaults to false." << std::endl;
	std::cout << "Set LAS2CSRS_DATA to point to ITRF DB and grid shift file." << std::endl;
}

int test(int argc, char **argv) {

	if(argc < 10) {
		std::cerr << "Too few arguments." << std::endl;
		return 1;
	}
	try {

		int i = 1;
		std::string ffrom(argv[++i]);
		double efrom = atof(argv[++i]);
		double eto = atof(argv[++i]);
		int fsrid = atoi(argv[++i]);
		int tsrid = atoi(argv[++i]);
		double x = atof(argv[++i]);
		double y = atof(argv[++i]);
		double z = atof(argv[++i]);

		Transformer trans(ffrom, efrom, eto, fsrid, tsrid);

		MemRaster<double> xx(1,1);
		MemRaster<double> yy(1,1);
		MemRaster<double> zz(1,1);

		xx.set(0, x);
		yy.set(0, y);
		zz.set(0, z);

		double bounds[6];

		trans.transformPoints(xx, yy, zz, 1, bounds);
		std::cout << std::setprecision(12) << xx[0] << " " << yy[0] << " " << zz[0] << std::endl;

	} catch(const char *e) {
		std::cerr << e << std::endl;
		return 1;
	}

	return 0;
}

int main(int argc, char **argv) {

	std::string comp("test");
	if(argc >= 2 && comp.compare(argv[1]) == 0) {
		return test(argc, argv);
	}

	if(argc < 8) {
		usage();
		return 1;
	}

	bool overwrite = false;

	try {

		int i = 1;
		for(; i < argc; ++i) {
			std::string arg(argv[i]);
			if(arg.c_str()[0] != '-') {
				break;
			} else if(arg == "-o") {
				overwrite = true;
			} else if(arg == "-v") {
				quiet = false;
			}
		}

		std::string srcfile(argv[i]);
		std::string dstdir(argv[++i]);
		std::string ffrom(argv[++i]);
		double efrom = atof(argv[++i]);
		double eto = atof(argv[++i]);
		int fsrid = atoi(argv[++i]);
		int tsrid = atoi(argv[++i]);

		if(!std::getenv(LAS2CSRS_DATA))
			setenv(LAS2CSRS_DATA, "..", 1);

		Transformer trans(ffrom, efrom, eto, fsrid, tsrid);

		trans.transformLas(srcfile, dstdir, overwrite);

	} catch(const std::exception &err) {
		std::cerr << err.what() << std::endl;
		usage();
		return 1;
	} catch(const std::string &err) {
		std::cerr << err << std::endl;
		usage();
		return 1;
	} catch(const char *err) {
		std::cerr << err << std::endl;
		usage();
		return 1;
	} catch(...) {
		std::cerr << "An unknown exception occurred." << std::endl;
		usage();
		return 1;
	}

	return 0;

}

