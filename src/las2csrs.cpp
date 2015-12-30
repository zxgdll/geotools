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

#include "Util.hpp"

#define PI 3.14159265358979323846
#define LAS2CSRS_DATA "LAS2CSRS_DATA"

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;
namespace las = liblas;


// Convert from arcsec to radians.
double _sec2rad(double x) {
	return x * 4.84813681 / 1000000000.0;
}

// Convert from radians to degrees.
double _deg(double x) {
	return x * 180.0 / PI;
}

double _rad(double x) {
	return x * PI / 180.0;
}

double _sq(double x) {
	return x*x;
}

double _binterp(float *grid, double c, double r, int c0, int r0, int c1, int r1, int width) {
	double x1 = (c1 - c) / (c1 - c0) * grid[r0 * width + c0] + (c - c0) / (c1 - c0) * grid[r0 * width + c1];
	double x2 = (c1 - c) / (c1 - c0) * grid[r1 * width + c0] + (c - c0) / (c1 - c0) * grid[r1 * width + c1];
	return (r1 - r) / (r1 - r0) * x1 + (r - r0) / (r1 - r0) * x2;
}

/**
 * Convert projected distance in mm to latlon in radians.
 * dx, dy  -- Distance in m.
 * lat     -- The latitude at which distances are computed
 * h       -- Ellipsoidal height.
 * a       -- Semi-major axis
 * e2      -- Eccentricity^2
 * count   -- Number of points.
 * dlat    -- The distance in rad (out).
 * dlon    -- The distance in rad (out).
 */
void _shift2latlon(double *dx, double *dy, double *lat, double *h, double a, double e2, int count, 
		double *dlat, double *dlon) {
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
	void interpolate(double *x, double *y, double *dx, double *dy, double *dz, int count) {
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
	Grid<float> xg;
	Grid<float> yg;
	Grid<float> zg;
	Grid<double> tg;
	int width;
	int height;
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
		this->efrom = efrom;
		this->eto = eto;
		this->fsrid = fsrid;
		this->tsrid = tsrid;
		this->ffrom.assign(ffrom);

		initProjections();
		loadHelmert();
		shiftGrid.load();
	}

	/**
	 * Transforms coordinate(s) from one available reference frame to NAD83(CSRS).
	 * x, y, z -- Coordinate arrays.
	 * count   -- The number of coordinates.
	 * bounds  -- The bounds of the transformed coordinate list.
	 */
	void transformPoints(Grid<double> &x, Grid<double> &y, Grid<double> &z, int count, double bounds[6]) {

		// Project to Cartesian 3D. (c)
		pj_transform(p1, p2, count, 1, x.grid(), y.grid(), z.grid());

		// Transform to csrs using Helmert (etc.) params. (c)
		epochTransform(x.grid(), y.grid(), z.grid(), count, this->efrom - 1997.0);

		// Only use the grid shift if the epoch changes.
		if(efrom != eto) {

			// Copy the coordinate arrays for transformation.
			Grid<double> x0(count, 1);
			Grid<double> y0(count, 1);
			Grid<double> z0(count, 1);
			memcpy(x0.grid(), x.grid(), sizeof(double) * count);
			memcpy(y0.grid(), y.grid(), sizeof(double) * count);
			memcpy(z0.grid(), z.grid(), sizeof(double) * count);

			// Initalize shift arrays.
			Grid<double> dx(count, 1);
			Grid<double> dy(count, 1);
			Grid<double> dz(count, 1);
			
			// Transform to latlon. (b)
			pj_transform(p2, p4, count, 1, x0.grid(), y0.grid(), z0.grid());
			
			// Interpolate shifts using latlon coords -- returns in mm. (d)
			shiftGrid.interpolate(x0.grid(), y0.grid(), dx.grid(), dy.grid(), dz.grid(), count);
			
			// Transform mm shifts to latlon
			Grid<double> dlat(count, 1);
			Grid<double> dlon(count, 1);

			// Get projection's spheroid props.
			double a, e2;
			pj_get_spheroid_defn(p4, &a, &e2);
			// Change the shift distance in projected coords to latlon.
			// This avoids the scale distortion associated with projection.
			_shift2latlon(dx.grid(), dy.grid(), y0.grid(), z0.grid(), a, e2, count, dlat.grid(), dlon.grid());

			double dt = this->eto - this->efrom;
			// Apply shifts to latlon coords.
			for(int i = 0; i < count; ++i) {
				*(x0.grid() + i) += *(dlon.grid() + i) * dt;
				*(y0.grid() + i) += *(dlat.grid() + i) * dt;
				*(z0.grid() + i) += *(dz.grid() + i) * dt;
			}
			
			// Transform latlon to target proj
			pj_transform(p4, p3, count, 1, x0.grid(), y0.grid(), z0.grid());
			
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
		
		const fs::path src(srcfile);
		const fs::path dst(dstdir);

		// Check the files for sanity.
		if(src == dst)
			throw std::string("Destination and source are the same: ") + dst.string() + std::string("==") + src.string();
		//if(!fs::is_directory(dst))
		//	throw std::string("Destination is a file: ") + dst.string();
		if(!fs::exists(dst)) {
			if(!fs::create_directory(dst))
				throw std::string("Failed to create output directory: ") + dst.string();
		}

		// Get the list of las files.
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
				if((overwrite || !fs::exists(dst / di->path().leaf())) && alg::ends_with(p, ext))
					files.push_back(di->path());
			}
		}
		if(files.size() == 0)
			throw "No matching files found.";

		std::cout << "Processing " << files.size() << " files." << std::endl;

		// Start
		las::WriterFactory wf;
		las::ReaderFactory rf;

		// The overall bounds: min x, max x, min y, max y, min z, max z
		double bounds[] = { FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX };

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

			Grid<double> x(count, 1);
			Grid<double> y(count, 1);
			Grid<double> z(count, 1);

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

	// Source reference frame.
	std::string ffrom;
	// From and to epochs.
	double efrom;
	double eto;
	// Projection objects.
	projPJ p1;
	projPJ p2;
	projPJ p3;
	projPJ p4;
	// SRIDs of the from and to CRSes.
	int fsrid;
	int tsrid;
	// The shift grid.
	ShiftGrid shiftGrid;
	// Transform parameters; loaded from the itrf file.
	double tx, ty, tz, dtx, dty, dtz;		// Shifts, rates.
	double rx, ry, rz, drx, dry, drz;		// Rotations, rates.
	double epoch;							// ITRF Transform epoch.
	double dt; 								// Time delta
	double d, dd; 							// Scale, scale rate.

	/**
		Transform the coordinate using the procedure listed in Craymer (2006).

		x, y, z -- 	The coordinate arrays.
		count 	-- 	The number of points.
	*/
	void epochTransform(double *x, double *y, double *z, int count, double dt) {
		double a0 = tx + dtx * dt;
		double a1 = ty + dty * dt;
		double a2 = tz + dtz * dt;
		double bsx = 1.0 + d + dd * dt;
		double b01 = -_sec2rad(rz + drz * dt);
		double b02 = _sec2rad(ry + dry * dt);
		double b10 = _sec2rad(rz + drz * dt);
		double b12 = -_sec2rad(rx + drx * dt);
		double b20 = -_sec2rad(ry + dry * dt);
		double b21 = _sec2rad(rx + drx * dt);
		for(;count;--count, ++x, ++y, ++z) {
			*x = a0 + bsx * *x + b01 * *y + b02 * *z;
			*y = a1 + b10 * *x + bsx * *y + b12 * *z;
			*z = a2 + b20 * *x + b21 * *y + bsx * *z;
		}
	}

	/**
	 * Initialize the proj4 projection objects.
	 */
	void initProjections() {
		// Initialize projections.
		if(!fsrid || !tsrid)
			throw "SRIDs are not set.";
		char str[128];
		sprintf(str, "+init=epsg:%u", fsrid);
		p1 = pj_init_plus(str);
		if(!p1) {
			const char *err = pj_strerrno(pj_errno);
			throw err;
		}
		sprintf(str, "+init=epsg:%u", tsrid);
		p3 = pj_init_plus(str);
		if(!p3){
			const char *err = pj_strerrno(pj_errno);
			throw err;
		}
		p2 = pj_init_plus("+init=epsg:4978");
		if(!p2){
			const char *err = pj_strerrno(pj_errno);
			throw err;
		}
		p4 = pj_init_plus("+init=epsg:4326");
		if(!p4){
			const char *err = pj_strerrno(pj_errno);
			throw err;
		}
	}

	/**
	 * Load the transformation database.
	 */
	void loadHelmert() {

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
			if(_ffrom == this->ffrom) {
				this->epoch = epoch;
				this->d = d / 1000000000.0; 		// Listed as ppb.
				this->dd = dd / 1000000000.0;
				this->tx = tx;
				this->ty = ty;
				this->tz = tz;
				this->rx = rx;
				this->ry = ry;
				this->rz = rz;
				this->dtx = dtx;
				this->dty = dty;
				this->dtz = dtz;
				this->drx = drx;
				this->dry = dry;
				this->drz = drz;
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
}

int test(int argc, char **argv) {

	if(argc < 10) {
		std::cerr << "Too few arguments." << std::endl;
		return 1;
	}

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

	Grid<double> xx(1,1);
	Grid<double> yy(1,1);
	Grid<double> zz(1,1);

	xx[0] = x;
	yy[0] = y;
	zz[0] = z;
	double bounds[6];

	trans.transformPoints(xx, yy, zz, 1, bounds);

	std::cout << std::setprecision(12) << xx[0] << " " << yy[0] << " " << zz[0] << std::endl;

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
		return 1;
	} catch(const std::string &err) {
		std::cerr << err << std::endl;
		return 1;
	} catch(...) {
		std::cerr << "An unknown exception occurred." << std::endl;
		return 1;
	}

	return 0;

}

