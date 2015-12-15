#!/usr/bin/env python

import pyproj
import os
import sys
import gdal
import subprocess
from math import sin, cos

PI = 3.14159265358979323846
FLT_MIN = -99999999.9
FLT_MAX = 99999999.9
ITRF_FILE = '../share/itrf.csv'
SHIFT_FILE = '../share/NAD83v6VG.tif'

# Convert from arcsec to radians.
def _sec2rad(x):
	return x * 4.84813681 / 1000000000.0

# Convert from radians to degrees.
def _deg(x):
	return x * 180.0 / PI

def _rad(x):
	return x * PI / 180.0

def _sq(x):
	return x*x

def _binterp(grid, c, r, c0, r0, c1, r1, width):
	'''
	// w -- x  c0 -- c1
	// |    |   |    |
	// y -- z  r0 -- r1
	'''
	a1 = ((c1 - c) * (r1 - r))
	g1 = grid[r0, c0]
	a2 = ((c - c0) * (r1 - r))
	g2 = grid[r0, c1]
	a3 = ((c1 - c) * (r - r0))
	g3 = grid[r1, c0]
	a4 = ((c - c0) * (r - r0))
	g4 = grid[r1, c1]
	return a1 * g1 + a2 * g2 + a3 * g3 + a4 * g4

def _shift2latlon(dx, dy, lat, h, a, e2):
	'''
	 * Convert projected distance in mm to latlon in radians.
	 * dx, dy  -- Distance in m.
	 * lat     -- The latitude at which distances are computed
	 * h       -- Ellipsoidal height.
	 * a       -- Semi-major axis
	 * e2      -- Eccentricity^2
	 * count   -- Number of points.
	 * dlat    -- The distance in rad (out).
	 * dlon    -- The distance in rad (out).
	'''
	dlat = [0.] * len(dx)
	dlon = [0.] * len(dx)
	for i in range(len(dx)):
		l = lat[i]
		m = a * (1.0 - e2) / (1.0 - e2 * _sq(sin(l))) ** (3.0/2.0) 	# Meridional radius of curvature.
		n = a / (1.0 - e2 * sin(l) ** 2) ** (1.0/2.0) 				# Parallel radius of curvature.
		r = n * cos(l) 												# Radius of parallel.
		dlon[i] = dx[i] / (r + h[i])
		dlat[i] = dy[i] / (m + h[i])
	return (dlon, dlat)

class ShiftGrid(object):
	'''
	 * Loads and interpolates the NAD83(CSRS) shift grid.
	'''

	def load(self):
		'''
		/**
		 * Load the shift grid raster and store the x, y, z shifts internally.
		 */
		 // TODO: Only load the portion of the grid necessary for the point cloud.
		 //       On the other hand, the grid is so small, who cares?
		'''
		ds = gdal.Open(SHIFT_FILE, gdal.GA_ReadOnly)
		if not ds:
			raise Exception("Failed to load shift grid.")
		xband = ds.GetRasterBand(1)
		yband = ds.GetRasterBand(2)
		zband = ds.GetRasterBand(3)
		self.width = xband.XSize
		self.height = xband.YSize
		self.xgrid = xband.ReadAsArray()
		self.ygrid = yband.ReadAsArray()
		self.zgrid = zband.ReadAsArray()
		self.transform = ds.GetGeoTransform()

	def interpolate(self, x, y, dx, dy, dz):
		'''
		/**
		 * Compute the the shifts in x, y and z at position x, y. x and y are radians,
		 * dx, dy and dz are m.
		 */
		'''
		transform = self.transform
		width = self.width
		height = self.height
		xgrid = self.xgrid
		ygrid = self.ygrid
		zgrid = self.zgrid
		for i in range(len(x)):
			c = ((_deg(x[i]) - transform[0]) / transform[1])
			r = ((_deg(y[i]) - transform[3]) / transform[5])
			c0 = int(c)
			r0 = int(r)
			c1 = c0 + 1
			r1 = r0 + 1
			if(c0 < 0):
				c0 = 0
			if(r0 < 0):
				r0 = 0
			if(c1 >= width):
				c1 = width - 1
			if(r1 >= height):
				r1 = height - 1
			dx[i] = _binterp(xgrid, c, r, c0, r0, c1, r1, width) / 1000.0
			dy[i] = _binterp(ygrid, c, r, c0, r0, c1, r1, width) / 1000.0
			dz[i] = _binterp(zgrid, c, r, c0, r0, c1, r1, width) / 1000.0


class Transformer(object):
	'''
	/**
	 * Performs the work of transforming coordinates from a given reference frame to NAD83(CSRS).
 	*/
 	'''

	def __init__(self, ffrom, efrom, eto, fsrid, tsrid):
		'''
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
	 	'''
		self.efrom = efrom
		self.eto = eto
		self.fsrid = fsrid
		self.tsrid = tsrid
		self.ffrom = ffrom.lower()

		self.initProjections()
		self.loadHelmert()
		self.shiftGrid = ShiftGrid()
		self.shiftGrid.load()


	def transformPoints(self, x, y, z):
		'''
		/**
		 * Transforms coordinate(s) from one available reference frame to NAD83(CSRS).
		 * a = UTM (m)
	 	 * c = a -> 4798 (geocentric, m) -> helmert
	 	 * b = c -> 4326 (geodetic, deg)
	 	 * d = b -> shifts (mm)
	 	 * e = d -> geocentric (sec)
	 	 * f = c + e
	 	 * g = f -> UTM (m)
	 	 *
		 * x, y, z -- Coordinate arrays.
		 * count   -- The number of coordinates.
		*/
		'''
		if len(x) != len(y) or len(x) != len(z):
			raise Exception('Point arrays must be the same length.')

		count = len(x)
		# Project to Cartesian 3D. (c)
		x, y, z = pyproj.transform(self.p1, self.p2, x, y, z, True)

		# Transform to csrs using Helmert (etc.) params. (c)
		self.epochTransform(x, y, z)

		# Only use the grid shift if the epoch changes.
		if self.efrom != self.eto:

			# Copy the coordinate arrays for transformation.
			x0 = x[:]
			y0 = y[:]
			z0 = z[:]

			# Initalize shift arrays.
			dx = [0.] * count
			dy = [0.] * count
			dz = [0.] * count
			
			# Transform to latlon. (b)
			x0, y0, z0 = pyproj.transform(self.p2, self.p4, x0, y0, z0, True)
			
			# Interpolate shifts using latlon coords -- returns in mm. (d)
			self.shiftGrid.interpolate(x0, y0, dx, dy, dz)
			
			# Transform mm shifts to latlon
			# Get projection's spheroid props.
			#double a, e2
			#pj_get_spheroid_defn(p4, &a, &e2)
			a, e2 = self.getEllipsProps("+init=EPSG:4326")

			# Change the shift distance in projected coords to latlon.
			# This avoids the scale distortion associated with projection.
			dlon, dlat = _shift2latlon(dx, dy, y0, z0, a, e2)

			dt = self.dt
			# Apply shifts to latlon coords.
			for i in range(count):
				x0[i] += dlon[i] * dt
				y0[i] += dlat[i] * dt
				z0[i] += dz[i] * dt
			
			# Transform latlon to target proj
			x, y, z = pyproj.transform(self.p4, self.p3, x0, y0, z0, True)

		else:

			# Reproject to the dest coordinates
			x, y, z = pyproj.transform(self.p2, self.p3, x, y, z, True)

		bounds = [FLT_MAX, FLT_MIN, FLT_MAX, FLT_MIN, FLT_MAX, FLT_MIN]
		# Expand the bounds for the new header.
		for i in range(count):
			if x[i] < bounds[0]: bounds[0] = x[i]
			if x[i] > bounds[1]: bounds[1] = x[i]
			if y[i] < bounds[2]: bounds[2] = y[i]
			if y[i] > bounds[3]: bounds[3] = y[i]
			if z[i] < bounds[4]: bounds[4] = z[i]
			if z[i] > bounds[5]: bounds[5] = z[i]

		#printf("xmin: %0.6f xmax: %0.6f ymin: %0.6f ymax: %0.6f zmin: %0.6f zmax: %0.6f\n",
		#		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5])

		return (x, y, z, bounds)

	def epochTransform(self, x, y, z):
		'''
		/**
		Transform the coordinate using the procedure listed in Craymer (2006).

		x, y, z -- 	The coordinate arrays.
		count 	-- 	The number of points.
		*/
		'''
		bsx = self.bsx
		b01 = self.b01
		b02 = self.b02
		b10 = self.b10
		b12 = self.b12
		b20 = self.b20
		b21 = self.b21
		a0 = self.a0
		a1 = self.a1
		a2 = self.a2
		for i in range(len(x)):
			x[i] = a0 + bsx * x[i] + b01 * y[i] + b02 * z[i]
			y[i] = a1 + b10 * x[i] + bsx * y[i] + b12 * z[i]
			z[i] = a2 + b20 * x[i] + b21 * y[i] + bsx * z[i]

	def getEllipsProps(self, proj):
		p = subprocess.Popen(['proj_ellips', proj], stdout=subprocess.PIPE)
		stdout, stderr = p.communicate()
		return map(float, stdout.split())

	def initProjections(self):
		# Initialize projections.
		if not self.fsrid or not self.tsrid:
			raise Exception("SRIDs are not set.")
		self.p1 = pyproj.Proj("+init=EPSG:%u" % (self.fsrid,))
		self.p3 = pyproj.Proj("+init=EPSG:%u" % (self.tsrid,))
		self.p2 = pyproj.Proj("+init=EPSG:4978")
		self.p4 = pyproj.Proj("+init=EPSG:4326")

	def loadHelmert(self):
		'''
		/**
		 * Load the transformation database.
		 */
		'''
		found = False
		with open(ITRF_FILE, "rU") as f:
			line = f.readline()
			while line:
				line = line.strip().split()
				if line[0] == self.ffrom:
					ffrom, fto = line[:2]
					epoch, tx, ty, tz, rx, ry, rz, d, dtx, dty, dtz, drx, dry, drz, dd = map(float, line[2:])
					dt = self.dt = self.eto - self.efrom
					self.d = d / 1000000000 		# Listed as ppb.
					self.dd = dd / 1000000000
					self.a0 = tx + dtx * dt
					self.a1 = ty + dty * dt
					self.a2 = tz + dtz * dt
					self.bsx = 1 + d + dd * dt
					self.b01 = -_sec2rad(rz + drz * dt)
					self.b02 = _sec2rad(ry + dry * dt)
					self.b10 = _sec2rad(rz + drz * dt)
					self.b12 = -_sec2rad(rx + drx * dt)
					self.b20 = -_sec2rad(ry + dry * dt)
					self.b21 = _sec2rad(rx + drx * dt)
					found = True
					break
				line = f.readline()

		if not found:
			raise Exception("Failed to find a transformation matching the parameters.")

if __name__ == '__main__':

	ffrom = sys.argv[1]
	efrom = float(sys.argv[2])
	eto = float(sys.argv[3])
	fsrid = int(sys.argv[4])
	tsrid = int(sys.argv[5])

	t = Transformer(ffrom, efrom, eto, fsrid, tsrid)
	x = []
	y = []
	z = []
	for i in range(6, len(sys.argv), 3):
		x.append(float(sys.argv[i]))
		y.append(float(sys.argv[i+1]))
		z.append(float(sys.argv[i+2]))

	x, y, z, bounds = t.transformPoints(x, y, z)
	print x, y, z, bounds
