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
		transform = (double *) malloc(6 * sizeof(double));
		ds = (GDALDataset *) GDALOpen(path, GA_ReadOnly);
		if(!ds)
			throw "Failed to load shift grid.";
		xband = ds->GetRasterBand(1);
		yband = ds->GetRasterBand(2);
		zband = ds->GetRasterBand(3);
		if(!xband || !yband || !zband) {
			GDALClose(ds);
			throw "Failed to retrieve shift bands.";
		}
		width = xband->GetXSize();
		height = yband->GetYSize();
		if(width <= 0 || height <= 0) {
			GDALClose(ds);
			throw "The dimensions of the shift grid are invalid.";
		}
		xgrid = (float *) malloc((unsigned long) (width * height * sizeof(float)));
		ygrid = (float *) malloc((unsigned long) (width * height * sizeof(float)));
		zgrid = (float *) malloc((unsigned long) (width * height * sizeof(float)));
		if(!xgrid || !ygrid || !zgrid) {
			free(xgrid);
			free(ygrid);
			free(zgrid);
			GDALClose(ds);
			throw "Failed to initialize shift grids.";
		}
		CPLErr xe = xband->RasterIO(GF_Read, 0, 0, width, height, xgrid, width, height, GDT_Float32, 0, 0);
		CPLErr ye = yband->RasterIO(GF_Read, 0, 0, width, height, ygrid, width, height, GDT_Float32, 0, 0);
		CPLErr ze = zband->RasterIO(GF_Read, 0, 0, width, height, zgrid, width, height, GDT_Float32, 0, 0);
		if(xe == CE_Failure || ye == CE_Failure || ze == CE_Failure) {
			free(xgrid);
			free(ygrid);
			free(zgrid);
			GDALClose(ds);
			throw "Failed to read shift grid.";
		}
		ds->GetGeoTransform(transform);
		GDALClose(ds);
	}

	/**
	 * Compute the the shifts in x, y and z at position x, y. x and y are radians,
	 * dx, dy and dz are m.
	 */
	void interpolate(double *x, double *y, double *dx, double *dy, double *dz, int count) {
		double c, r;
		int c0, c1, r0, r1;
		for(;count;--count, ++x, ++y, ++dx, ++dy, ++dz) {
			c = ((_deg(*x) - transform[0]) / transform[1]);
			r = ((_deg(*y) - transform[3]) / transform[5]);
			c0 = (int) c;
			r0 = (int) r;
			c1 = c0 + 1;
			r1 = r0 + 1;
			if(c0 < 0) c0 = 0;
			if(r0 < 0) r0 = 0;
			if(c1 >= width) c1 = width-1;
			if(r1 >= height) r1 = height-1;
			*dx = _binterp(xgrid, c, r, c0, r0, c1, r1, width) / 1000.0;
			*dy = _binterp(ygrid, c, r, c0, r0, c1, r1, width) / 1000.0;
			*dz = _binterp(zgrid, c, r, c0, r0, c1, r1, width) / 1000.0;
		}
	}

	~ShiftGrid() {
		free(xgrid);
		free(ygrid);
		free(zgrid);
		free(transform);
	}

private:

	float *xgrid;
	float *ygrid;
	float *zgrid;
	double *transform;
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
	void transformPoints(double *x, double *y, double *z, int count, double bounds[6]) {

		// Project to Cartesian 3D. (c)
		pj_transform(p1, p2, count, 1, x, y, z);

		// Transform to csrs using Helmert (etc.) params. (c)
		epochTransform(x, y, z, count, this->efrom - 1997.0);

		// Only use the grid shift if the epoch changes.
		if(efrom != eto) {

			// Copy the coordinate arrays for transformation.
			double *x0 = (double *) malloc(sizeof(double) * count);
			double *y0 = (double *) malloc(sizeof(double) * count);
			double *z0 = (double *) malloc(sizeof(double) * count);
			memcpy(x0, x, sizeof(double) * count);
			memcpy(y0, y, sizeof(double) * count);
			memcpy(z0, z, sizeof(double) * count);

			// Initalize shift arrays.
			double *dx = (double *) malloc(sizeof(double) * count);
			double *dy = (double *) malloc(sizeof(double) * count);
			double *dz = (double *) malloc(sizeof(double) * count);
			
			// Transform to latlon. (b)
			pj_transform(p2, p4, count, 1, x0, y0, z0);
			
			// Interpolate shifts using latlon coords -- returns in mm. (d)
			shiftGrid.interpolate(x0, y0, dx, dy, dz, count);
			
			// Transform mm shifts to latlon
			double *dlat = (double *) malloc(sizeof(double) * count);
			double *dlon = (double *) malloc(sizeof(double) * count);
			// Get projection's spheroid props.
			double a, e2;
			pj_get_spheroid_defn(p4, &a, &e2);
			// Change the shift distance in projected coords to latlon.
			// This avoids the scale distortion associated with projection.
			_shift2latlon(dx, dy, y0, z0, a, e2, count, dlat, dlon);
			free(dx);
			free(dy);

			double dt = this->eto - this->efrom;
			// Apply shifts to latlon coords.
			for(int i = 0; i < count; ++i) {
				*(x0 + i) += *(dlon + i) * dt;
				*(y0 + i) += *(dlat + i) * dt;
				*(z0 + i) += *(dz + i) * dt;
			}
			free(dlat);
			free(dlon);
			free(dz);
			
			// Transform latlon to target proj
			pj_transform(p4, p3, count, 1, x0, y0, z0);
			
			// Assign the shifted coords to the output arrays.
			memcpy(x, x0, sizeof(double) * count);
			memcpy(y, y0, sizeof(double) * count);
			memcpy(z, z0, sizeof(double) * count);
			free(x0);
			free(y0);
			free(z0);

		} else {

			// Reproject to the dest coordinates
			pj_transform(p2, p3, count, 1, x, y, z);

		}

		// Expand the bounds for the new header.
		for(;count;--count, ++x, ++y, ++z) {
			if(*x < bounds[0]) bounds[0] = *x;
			if(*x > bounds[1]) bounds[1] = *x;
			if(*y < bounds[2]) bounds[2] = *y;
			if(*y > bounds[3]) bounds[3] = *y;
			if(*z < bounds[4]) bounds[4] = *z;
			if(*z > bounds[5]) bounds[5] = *z;
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
		// TODO: Ignore completed files.
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

			double *x = (double *) malloc(sizeof(double) * count);
			double *y = (double *) malloc(sizeof(double) * count);
			double *z = (double *) malloc(sizeof(double) * count);

			// Iterate over the points.
			for(int i = 0; r.ReadNextPoint(); ++i) {
				las::Point pt = r.GetPoint();
				*(x + i) = pt.GetX();
				*(y + i) = pt.GetY();
				*(z + i) = pt.GetZ();
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
				pt.SetX(*(x + i));
				pt.SetY(*(y + i));
				pt.SetZ(*(z + i));
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
		sprintf(str, "+init=EPSG:%u", fsrid);
		p1 = pj_init_plus(str);
		if(!p1) {
			const char *err = pj_strerrno(pj_errno);
			throw err;
		}
		sprintf(str, "+init=EPSG:%u", tsrid);
		p3 = pj_init_plus(str);
		if(!p3){
			const char *err = pj_strerrno(pj_errno);
			throw err;
		}
		p2 = pj_init_plus("+init=EPSG:4978");
		if(!p2){
			const char *err = pj_strerrno(pj_errno);
			throw err;
		}
		p4 = pj_init_plus("+init=EPSG:4326");
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

	double xx[1] = {x};
	double yy[1] = {y};
	double zz[1] = {z};
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

		int i = 0;
		for(; i < argc; ++i) {
			std::string arg(argv[i]);
			if(arg.c_str()[0] != '-') {
				break;
			} else if(arg == "-o") {
				overwrite = true;
			}
		}

		std::string srcfile(argv[++i]);
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

	} catch(char const *err) {
		std::cerr << err << std::endl;
		return 1;
	} catch(const std::string &err) {
		std::cerr << err << std::endl;
	}

	return 0;

}

