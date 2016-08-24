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
#include "util.hpp"
#include "raster.hpp"

#define LAS2CSRS_DATA "LAS2CSRS_DATA"


#include "raster.cpp"

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;
namespace las = liblas;

using namespace geotools::util;
using namespace geotools::raster;

namespace geotools {

	namespace geodesy {

		namespace util {

			// Binary interpolation.
			static double _binterp(float *grid, double c, double r, int c0, int r0, int c1, int r1, int width) {
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
			static void _shift2latlon(Grid<double> &gdx, Grid<double> &gdy, Grid<double> &glat, Grid<double> &gh,
				double a, double e2, int count, Grid<double> &gdlat, Grid<double> &gdlon) {

				double *dx = gdx.grid();
				double *dy = gdy.grid();
				double *lat = glat.grid();
				double *dlat = gdlat.grid();
				double *dlon = gdlon.grid();
				double *h = gh.grid();

				for(int i = 0; i < count; ++i) {
					double l = *(lat + i);
					double m = a * (1.0 - e2) / pow((1.0 - e2 * g_sq(sin(l))), 3.0/2.0); 	// Meridional radius of curvature.
					double n = a / pow((1.0 - e2 * g_sq(sin(l))), 1.0/2.0); 					// Parallel radius of curvature.
					double r = n * cos(l); 													// Radius of parallel.
					*(dlon + i) = *(dx + i) / (r + *(h + i));
					*(dlat + i) = *(dy + i) / (m + *(h + i));		
				}
			}

			// Miliarcseconds to radians.
			inline double _mas2rad(double x) {
				return (x * 4.84813681 / 1000000000.0);
			}


			/**
			 * Loads and interpolates the NAD83(CSRS) shift grid.
			 */
			class ShiftGrid {
			private:
				MemRaster<float> xg;
				MemRaster<float> yg;
				MemRaster<float> zg;
				double tg[6];
				int width;
				int height;

			public:

				/**
				 * Load the shift grid raster and store the x, y, z shifts internally.
				 */
				void load() {

					GDALAllRegister();
					GDALDataset *ds;
					GDALRasterBand *xband, *yband, *zband;

					// Try to load from the adjacent share path
					ds = (GDALDataset *) GDALOpen("../share/NAD83v6VG.tif", GA_ReadOnly);
					// Otherwise try to load from the configured path.
					if(!ds) {
						char path[PATH_MAX];
						sprintf(path, "%s%s", std::getenv("LAS2CSRS_DATA"), "/NAD83v6VG.tif");
						ds = (GDALDataset *) GDALOpen(path, GA_ReadOnly);
					}
					if(!ds)
						g_runerr("Failed to load shift grid. Set the LAS2CSRS_DATA variable to point to the itrf database and grid shift file.");

					xband = ds->GetRasterBand(1);
					yband = ds->GetRasterBand(2);
					zband = ds->GetRasterBand(3);
					if(!xband || !yband || !zband)
						g_runerr("Failed to retrieve shift bands.");

					width = xband->GetXSize();
					height = yband->GetYSize();
					if(width <= 0 || height <= 0)
						g_runerr("The dimensions of the shift grid are invalid.");

					xg.init(width, height);
					yg.init(width, height);
					zg.init(width, height);

					CPLErr xe = xband->RasterIO(GF_Read, 0, 0, width, height, xg.grid(), width, height, GDT_Float32, 0, 0);
					CPLErr ye = yband->RasterIO(GF_Read, 0, 0, width, height, yg.grid(), width, height, GDT_Float32, 0, 0);
					CPLErr ze = zband->RasterIO(GF_Read, 0, 0, width, height, zg.grid(), width, height, GDT_Float32, 0, 0);

					if(xe == CE_Failure || ye == CE_Failure || ze == CE_Failure)
						g_runerr("Failed to read shift grid.");

					ds->GetGeoTransform(tg);
				}

				/**
				 * Compute the the shifts in x, y and z at position x, y. x and y are radians,
				 * dx, dy and dz are m.
				 */
				void interpolate(Grid<double> &gx, Grid<double> &gy, Grid<double> &gdx, 
					Grid<double> &gdy, Grid<double> &gdz) {

					int count = gx.size();
					double *x = gx.grid();
					double *y = gy.grid();
					double *dx = gdx.grid();
					double *dy = gdy.grid();
					double *dz = gdz.grid();

					double c, r;
					int c0, c1, r0, r1;
					for(;count;--count, ++x, ++y, ++dx, ++dy, ++dz) {
						c = ((g_deg(*x) - tg[0]) / tg[1]);
						r = ((g_deg(*y) - tg[3]) / tg[5]);
						c0 = (int) c;
						r0 = (int) r;
						c1 = c0 + 1;
						r1 = r0 + 1;
						if(c0 < 0) c0 = 0;
						if(r0 < 0) r0 = 0;
						if(c1 >= width) c1 = width-1;
						if(r1 >= height) r1 = height-1;
						*dx = _binterp(xg.grid(), c, r, c0, r0, c1, r1, width) / 1000.0; // Convert to m
						*dy = _binterp(yg.grid(), c, r, c0, r0, c1, r1, width) / 1000.0;
						*dz = _binterp(zg.grid(), c, r, c0, r0, c1, r1, width) / 1000.0;
					}
				}

				~ShiftGrid() {
				}
			};

			class Params {
			public:
				// Source reference frame.
				std::string ffrom;
				// From and to epochs.
				double efrom;
				double eto;
				// To and from SRSes.
				las::SpatialReference *fromSRS;
				las::SpatialReference *toSRS;
				// Transform parameters; loaded from the itrf file.
				double tx, ty, tz, dtx, dty, dtz;		// Shifts, rates.
				double rx, ry, rz, drx, dry, drz;		// Rotations, rates.
				double epoch;							// ITRF Transform epoch.
				//double dt; 							// Time delta
				double d, dd; 							// Scale, scale rate.

				Params() : 
					ffrom(""),
					efrom(0),
					eto(0),
					fromSRS(nullptr), 
					toSRS(nullptr),
					tx(0), ty(0), tz(0), dtx(0), dty(0), dtz(0),
					rx(0), ry(0), rz(0), drx(0), dry(0), drz(0),
					epoch(0), 
					d(0), dd(0) {}


				~Params() {
					if(fromSRS)
						delete fromSRS;
					if(toSRS)
						delete toSRS;
				}
			};

		} // Util

		using namespace geotools::geodesy::util;

		/**
		 * Performs the work of transforming coordinates from a given reference frame to NAD83(CSRS).
		 */
		class Transformer {
		private:
			// Projection objects.
			projPJ projFrom;
			projPJ projECEF;
			projPJ projTo;
			projPJ projGeog;
			// The shift grid.
			ShiftGrid shiftGrid;
			// Transformation params;
			Params params;     

			/**
				Transform the coordinate using the procedure listed in Craymer (2006).

				x, y, z -- 	The coordinate arrays.
				count 	-- 	The number of points.
			*/
			void epochTransform(Params &params, Grid<double> &gx, Grid<double> &gy, Grid<double> &gz, int count, double dt) {

				g_trace("epochTransform");

				double *x = gx.grid();
				double *y = gy.grid();
				double *z = gz.grid();

				double Txt = params.tx + params.dtx * dt;            // Translation in X plus X velocity * time
				double Tyt = params.ty + params.dty * dt;
				double Tzt = params.tz + params.dtz * dt;
				double DSt = params.d + params.dd * dt;         	// Scale plus delta scale * time
				double Rxt = _mas2rad(params.rx + params.drx * dt); // Rotation in X plus velocity * time; seconds to radians.
				double Ryt = _mas2rad(params.ry + params.dry * dt);
				double Rzt = _mas2rad(params.rz + params.drz * dt);

				DSt += 1.0; // Scale is w/r/t 1.

				g_trace("  > Txt: " << Txt << "; Tyt: " << Tyt << "; Tzt: " << Tzt << "\n" \
					<< "  > Rxt: " << Rxt << "; Ryt: " << Ryt << "; Rzt: " << Rzt << "\n" \
					<< "  > DSt: " << DSt << "\n" \
					<< "  > drx: " << params.drx << "; dry: " << params.dry << "; drz: " << params.drz << "\n");

				for(;count;--count, ++x, ++y, ++z) {
					*x = Txt + (DSt * *x) + (-Rzt * *y) + (Ryt * *z);
					*y = Tyt + (Rzt * *x) + (DSt * *y) + (-Rxt * *z);
					*z = Tzt + (-Ryt * *x) + (Rxt * *y) + (DSt * *z);
				}
			}

			/**
			 * Initialize the proj4 projection objects.
			 */
			void initProjections() {

				g_trace("initProjection");

				// Initialize projections.
				if(!params.fromSRS || !params.toSRS)
					g_argerr("SRSes are not set.");
				g_trace("From projection: " << params.fromSRS->GetProj4());
				projFrom = pj_init_plus(params.fromSRS->GetProj4().c_str());
				if(!projFrom)
					g_argerr(pj_strerrno(pj_errno));
				g_trace("To projection: " << params.toSRS->GetProj4());
				projTo = pj_init_plus(params.toSRS->GetProj4().c_str());
				if(!projTo)
					g_argerr(pj_strerrno(pj_errno));
				projECEF = pj_init_plus("+proj=geocent +ellps=GRS80 +units=m +no_defs");
				if(!projECEF)
					g_argerr(pj_strerrno(pj_errno));
				projGeog = pj_init_plus("+proj=latlon +ellps=GRS80 +no_defs");
				if(!projGeog)
					g_argerr(pj_strerrno(pj_errno));
			}

			/**
			 * Load the transformation database.
			 */
			void loadHelmert(Params &p) {

				g_trace("loadHelmert " << p.ffrom);

				std::string ffrom;
				std::string fto;
				float epoch, tx, ty, tz, rx, ry, rz, d, dtx, dty, dtz, drx, dry, drz, dd;
				bool found = false;
				std::ifstream f;
				f.open("../share/itrf.csv");
				if(!f.is_open()) {
					char path[PATH_MAX];
					sprintf(path, "%s%s", std::getenv(LAS2CSRS_DATA), "/itrf.csv");
					f.open(path);
				}
				if(!f.is_open()) 
					g_runerr("Failed to open itrf database. If necessary, set LAS2CSRS_DATA environment variable.");
				std::string line;
				while(std::getline(f, line)) {
					if(line[0] == '/' || line[0] == ' ')
						continue;
					std::stringstream ls(line);
					ls >> ffrom >> fto >> epoch >> tx >> ty >> tz >> rx >> ry >> rz >> d 
						>> dtx >> dty >> dtz >> drx >> dry >> drz >> dd;

					g_trace(" -- Checking: " << ffrom);
					if(ffrom == p.ffrom) {
						g_trace(" -- Found entry for " << ffrom);
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

				if(!found)
					g_argerr("Failed to find a transformation matching the parameters.");

				g_trace(" -- Params: dtx: " << p.dtx << "; dty: " << p.dty << "; dtz: " << p.dtz);
				g_trace("            drx: " << p.drx << "; dry: " << p.dry << "; drz: " << p.drz);

			}

		public:

			/**
			 * Prepares the Helmert transformation parameters for the given transformation.
			 * These will be used in later method calls.
			 * The constructor will load an process the grid shift file and the transformation database.
			 * ffrom  -- The name of the reference frame, e.g. 'itrf90'
			 * efrom  -- The epoch of data collection (decimal years), e.g. 1994.2
			 * eto    -- The target epoch. One might select 1997.0 for BC or 2002.0 for Albera (etc.)
			 * fsrid  -- The SRID (EPSG code) of the source.
			 * tsrid  -- The SRID of the destination. This is the code for the UTM zone under
			 *           NAD83(CSRS), e.g. 2956 for UTM12N.
			 * fvsrid -- The input geoid model. This presumes input orthometric heights and the 
			 *           existence of a geoid model in the proj database. If the input is ellipsoidal, leave
			 *           as empty. (Example input: HT2_0.gtx for CGVD28)
			 * tvsrid -- The output geoid model. If ellipsoidal output is desired, leave as empty.
			 */
			Transformer(const std::string &ffrom, float efrom, float eto, const std::string &fromSRS, const std::string &toSRS) {

				if(ffrom.empty())
					g_argerr("A source reference frame is required.");
				if(efrom < 1980)
					g_argerr("A source epoch is required (should be 1980 or higher.)");
				if(eto < 1980)
					g_argerr("A destination epoch is required (should be 1980 or higher.)");
				if(toSRS.empty())
					g_argerr("A destination SRS is required. Include geoid information if orthometric output is required.");
				if(fromSRS.empty())
					g_warn("Attempting to load source SRS information from LAS files.");

				g_trace("Transformer:\n -- src ref frame: " << ffrom << "; src epoch: " << efrom \
					<< "; dst epoch: " << eto << " -- from srs: " << fromSRS << "; to srs: ");

				params.efrom = efrom;
				params.eto = eto;
				params.ffrom.assign(ffrom);
				if(!fromSRS.empty()) {
					params.fromSRS = new las::SpatialReference();
					params.fromSRS->SetFromUserInput(fromSRS);
				}
				if(!toSRS.empty()) {
					params.toSRS = new las::SpatialReference();
					params.toSRS->SetFromUserInput(toSRS);
				}

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

				// 1) Transform original points to their ECEF representation (Cartesian) at the original epoch.
				// 2) Deduce the transformation parameters to NAD83 @ 1997.0.
				// 3) Transform from NAD83 forward to target epoch.

				g_trace("transformPoints" \
					<< " -- Original: " << x.grid()[0] << ", " << y.grid()[0] << ", " << z.grid()[0]);

				// 1) Project to ECEF in the original reference frame.
				pj_transform(projFrom, projECEF, count, 1, x.grid(), y.grid(), z.grid());
				g_trace(" -- ECEF (Original): " << x.grid()[0] << ", " << y.grid()[0] << ", " << z.grid()[0]);

				// 2) Transform to NAD83 @ 1997.
				epochTransform(params, x, y, z, count, params.efrom - params.epoch);
				g_trace(" -- ECEF (CSRS): " << x.grid()[0] << ", " << y.grid()[0] << ", " << z.grid()[0]);

				// Only use the grid shift if the epoch changes.
				if(params.efrom != params.eto) {

					// 3) Transform from NAD83 @ 1997 to NAD83 at target epoch.

					// Copy the coordinate arrays for transformation.
					MemRaster<double> x0(count, 1);
					MemRaster<double> y0(count, 1);
					MemRaster<double> z0(count, 1);
					memcpy(x0.grid(), x.grid(), sizeof(double) * count);
					memcpy(y0.grid(), y.grid(), sizeof(double) * count);
					memcpy(z0.grid(), z.grid(), sizeof(double) * count);

					// Transform from Cartesian (CSRS) to latlon.
					pj_transform(projECEF, projGeog, count, 1, x0.grid(), y0.grid(), z0.grid());

					g_trace(" -- Lat Lon (CSRS): " << g_deg(x0.grid()[0]) << ", " << g_deg(y0.grid()[0]) << ", " << z0.grid()[0]);

					// Initalize shift arrays.
					MemRaster<double> dx(count, 1);
					MemRaster<double> dy(count, 1);
					MemRaster<double> dz(count, 1);
					
					// Interpolate shifts using latlon coords -- returns in mm. (d)
					shiftGrid.interpolate(x0, y0, dx, dy, dz);
					
					g_trace(" -- Grid Velocities: " << dx.grid()[0] << ", " << dy.grid()[0] << ", " << dz.grid()[0]);

					// Transform mm shifts to latlon
					MemRaster<double> dlat(count, 1);
					MemRaster<double> dlon(count, 1);

					// Get projection's spheroid props.
					double a, e2;
					pj_get_spheroid_defn(projTo, &a, &e2);

					// Get angular shifts from grid shift (mm).
					_shift2latlon(dx, dy, y0, z0, a, e2, count, dlat, dlon);

					g_trace(" -- Lat lon shifts: " << dlat.grid()[0] << ", " << dlon.grid()[0]);

					// Good to here...
					
					double dt = params.eto - params.efrom;
					// Apply shifts to latlon coords.
					for(int i = 0; i < count; ++i) {
						*(x0.grid() + i) += *(dlon.grid() + i) * dt;
						*(y0.grid() + i) += *(dlat.grid() + i) * dt;
						*(z0.grid() + i) += *(dz.grid() + i) * dt;
					}
					
					g_trace(" -- Lat lon shifted: " << g_deg(x0.grid()[0]) << ", " << g_deg(y0.grid()[0]) << ", " << z0.grid()[0]);

					// Transform latlon to target proj
					pj_transform(projGeog, projTo, count, 1, x0.grid(), y0.grid(), z0.grid());
					
					// Assign the shifted coords to the output arrays.
					memcpy(x.grid(), x0.grid(), sizeof(double) * count);
					memcpy(y.grid(), y0.grid(), sizeof(double) * count);
					memcpy(z.grid(), z0.grid(), sizeof(double) * count);

				} else {

					// Reproject to the dest coordinates
					pj_transform(projECEF, projTo, count, 1, x.grid(), y.grid(), z.grid());

				}

				g_trace(" -- Final: " << x.grid()[0] << ", " << y.grid()[0] << ", " << z.grid()[0]);

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
			int transformLas(std::vector<std::string> &srcfiles, std::string &dstdir, bool overwrite) {
				
				g_trace("transformLas");

				if(!srcfiles.size())
					g_argerr("No source files provided.");

				fs::path dst(dstdir);
				if(!fs::exists(dst)) {
					if(!fs::create_directory(dst))
						g_argerr("Failed to create output directory: " << dst.string());
				}

				g_trace("Processing " << srcfiles.size() << " files.");

				// Start
				las::WriterFactory wf;
				las::ReaderFactory rf;

				// The overall bounds: min x, max x, min y, max y, min z, max z
				double bounds[] = { G_DBL_MAX_POS, G_DBL_MAX_NEG, G_DBL_MAX_POS, G_DBL_MAX_NEG, G_DBL_MAX_POS, G_DBL_MAX_NEG };

				for(std::string filename:srcfiles) {

					g_trace("Processing file " << filename);

					// Open the source file.
					std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
					las::Reader r = rf.CreateWithStream(in);
					las::Header h = r.GetHeader();

					// If no spatial reference has been set for the input, try
					// to get it from the first las file.
					if(!params.fromSRS) {
						params.fromSRS = new las::SpatialReference(h.GetSRS());
						if(!params.fromSRS)
							g_argerr("SRS not provided and not found on first LAS file.");
					}
					
					// Check the files for sanity.
					const fs::path srcfile(filename);
					const fs::path dstfile(dst / srcfile.leaf());
					if(srcfile == dstfile)
						g_argerr("Destination and source are the same: " << dstfile.string() << "==" << srcfile.string());

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
					std::ofstream out(dstfile.c_str(), std::ios::out | std::ios::binary);
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
				pj_free(projFrom);
				pj_free(projECEF);
				pj_free(projTo);
				pj_free(projGeog);
			}

		};

	} // geodesy

} // geotools


void usage() {
	std::cerr   << "Usage: las2csrs <options> <las files>\n\n"
			    << "This program converts coordinates from LAS files from any reference frame to NAD83(CSRS),\n"
			    << "between any two epochs.\n\n"
			    << "If orthometric heights are used, be sure to provide SRSes with geoid parameters, or that the source\n"
			    << "files contain such information. SRSes are entered in the form, 'epsg:<horizontal code>+<vertical code>'.\n"
			    << "The + is only required if orthometric heights are desired.\n\n"

			    << " -o     Overwrite existing files. Defaults to false.\n"
				<< " -v     Verbose output.\n"
				<< " -d 	Destination folder. Required.\n"
				<< " -fs    The source SRS. If left out, will use the first las file's SRS. If the file has no SRS, will fail.\n"
				<< " -ts    The destination SRS. Required.\n"
				<< " -fe    The source epoch. Required.\n"
				<< " -te    The destination epoch. Required.\n"
				<< " -f     The source reference frame. Required.\n\n"

			    << "Set LAS2CSRS_DATA to point to ITRF DB and grid shift file.\n";
}

int main(int argc, char **argv) {

	try {

		bool overwrite = false;
		std::vector<std::string> srcfiles;
		std::string dstdir;
		std::string ffrom;
		double efrom;
		double eto;
		std::string toSRS;
		std::string fromSRS;

		for(int i = 1; i < argc; ++i) {
			std::string arg(argv[i]);
			if(arg == "-o") {
				overwrite = true;
			} else if(arg == "-v") {
				g_loglevel(1);
			} else if(arg == "-d") {
				dstdir.assign(argv[++i]);
			} else if(arg == "-fs") {
				fromSRS.assign(argv[++i]);
			} else if(arg == "-ts") {
				toSRS.assign(argv[++i]);
			} else if(arg == "-fe") {
				efrom = atof(argv[++i]);
			} else if(arg == "-te") {
				eto = atof(argv[++i]);
			} else if(arg == "-f") {
				ffrom.assign(argv[++i]);
			} else {
				srcfiles.push_back(argv[i]);
			}
		}

		if(!std::getenv(LAS2CSRS_DATA))
			setenv(LAS2CSRS_DATA, "..", 1);

		geotools::geodesy::Transformer trans(ffrom, efrom, eto, fromSRS, toSRS);
		trans.transformLas(srcfiles, dstdir, overwrite);

	} catch(const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;

}

