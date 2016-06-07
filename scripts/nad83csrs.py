#!/usr/bin/env python

import pyproj
import os
import sys
import gdal
import subprocess
from math import sin, cos
from numpy import dtype

PI = 3.14159265358979323846
FLT_MIN = -99999999.9
FLT_MAX = 99999999.9
ITRF_FILE = '../share/itrf.csv'
SHIFT_FILE = '../share/NAD83v6VG.tif'

# Convert from arcsec to radians.
def _sec2rad(x):
	return x * 4.84813681 / 1000000000.

# Convert from radians to degrees.
def _deg(x):
	return x * 180. / PI

def _rad(x):
	return x * PI / 180.

# TODO: This would probably be more accurate (i.e. closer to TRX) with bicubic?
def _binterp(grid, c, r, c0, r0, c1, r1, width):
	'''
	Bilinear interpolation of shift grid.
	'''
	x1 = (c1 - c) / (c1 - c0) * grid[r0, c0] + (c - c0) / (c1 - c0) * grid[r0, c1]
	x2 = (c1 - c) / (c1 - c0) * grid[r1, c0] + (c - c0) / (c1 - c0) * grid[r1, c1]
	return (r1 - r) / (r1 - r0) * x1 + (r - r0) / (r1 - r0) * x2

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
		m = a * (1. - e2) / (1. - e2 * sin(l) ** 2.) ** (3. / 2.) 	# Meridional radius of curvature.
		n = a / (1. - e2 * sin(l) ** 2.) ** (1. / 2.) 				# Parallel radius of curvature.
		r = n * cos(l) 												# Radius of parallel.
		dlon[i] = dx[i] / (r + h[i])
		dlat[i] = dy[i] / (m + h[i])
	return (dlon, dlat)

class ShiftGrid(object):
	'''
	Loads and interpolates the NAD83(CSRS) shift grid.
	'''

	def load(self):
		'''
		Load the shift grid raster and store the x, y, z shifts internally.
		TODO: Only load the portion of the grid necessary... On the other hand, 
		the grid is so small, who cares?
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
		Compute the the shifts in x, y and z at position x, y. x and y are radians,
		dx, dy and dz are m.
		'''
		transform = self.transform
		width = self.width
		height = self.height
		xgrid = self.xgrid
		ygrid = self.ygrid
		zgrid = self.zgrid
		for i in range(len(x)):
			_dx = _deg(x[i])
			_dy = _deg(y[i])
			c = ((_dx - transform[0]) / transform[1])
			r = ((_dy - transform[3]) / transform[5])
			c0 = int(c)
			r0 = int(r)
			c1 = c0 + 1
			r1 = r0 + 1
			if(c0 < 0): c0 = 0
			if(r0 < 0): r0 = 0
			if(c1 >= width): c1 = width - 1
			if(r1 >= height): r1 = height - 1
			dx[i] = _binterp(xgrid, c, r, c0, r0, c1, r1, width) / 1000.
			dy[i] = _binterp(ygrid, c, r, c0, r0, c1, r1, width) / 1000.
			dz[i] = _binterp(zgrid, c, r, c0, r0, c1, r1, width) / 1000.


class Params(object):

	def __init__(self):
		self.ffrom = None
		self.efrom = None
		self.eto = None
		self.fromSRS = None
		self.toSRS = None
		self.tx = None
		self.ty = None
		self.tz = None
		self.rx = None
		self.ry = None
		self.rz = None
		self.drx = None
		self.dry = None
		self.drz = None
		self.dtx = None
		self.dty = None
		self.dtz = None
		self.epoch = None
		self.d = None
		self.dd = None
		self.epoch = None

class Transformer(object):
	'''
	Performs the work of transforming coordinates from a given reference frame to NAD83(CSRS).
 	'''

	def __init__(self, ffrom, efrom, eto, fromSRS, toSRS):
		'''
	 	Prepares the Helmert transformation parameters for the given transformation.
		These will be used in later method calls.
		The constructor will load an process the grid shift file and the transformation database.
		ffrom -- The name of the reference frame, e.g. 'itrf90'
		efrom -- The epoch of data collection (decimal years), e.g. 1994.2
		eto   -- The target epoch. One might select 1997.0 for BC or 2002.0 for Albera (etc.)
		fromSRS -- The proj string for the source.
		toSRS -- The proj string of the destination. This is the code for the UTM zone under
		         NAD83(CSRS), e.g. 2956 for UTM12N.
	 	'''

	 	self.a = None

	 	self.params = Params()
		self.params.efrom = efrom
		self.params.eto = eto
		self.params.fromSRS = fromSRS
		self.params.toSRS = toSRS
		self.ffrom = ffrom.lower()

		self.initProjections()
		self.loadHelmert()
		self.shiftGrid = ShiftGrid()
		self.shiftGrid.load()


	def transformPoints(self, x, y, z):
		'''
		Transforms coordinate(s) from one available reference frame to NAD83(CSRS).
	 	x, y, z -- Coordinate arrays.
		'''
		if len(x) != len(y) or len(x) != len(z):
			raise Exception('Point arrays must be the same length.')

		count = len(x)
		# Project to Cartesian 3D. (c)
		x, y, z = pyproj.transform(self.projFrom, self.projECEF, x, y, z, True)
	
		# Transform to csrs using Helmert (etc.) params. (c)
		self.epochTransform(x, y, z, self.params.efrom - self.params.epoch)
	
		# Only use the grid shift if the epoch changes.
		if self.params.efrom != self.params.eto:

			# Copy the coordinate arrays for transformation.
			x0 = x[:]
			y0 = y[:]
			z0 = z[:]

			# Initalize shift arrays.
			dx = [0.] * count
			dy = [0.] * count
			dz = [0.] * count
			
			# Transform to latlon. (b)
			x0, y0, z0 = pyproj.transform(self.projECEF, self.projGeog, x0, y0, z0, True)
			
			# Interpolate shifts using latlon coords -- returns in mm. (d)
			self.shiftGrid.interpolate(x0, y0, dx, dy, dz)

			# Transform mm shifts to latlon
			# Get GRS80 spheroid props.
			a, e2 = self.getEllipsProps(self.projTo)

			# Change the shift distance in projected coords to latlon.
			# This avoids the scale distortion associated with projection.
			dlon, dlat = _shift2latlon(dx, dy, y0, z0, a, e2)

			dt = self.params.eto - self.params.efrom
			# Apply shifts to latlon coords.
			for i in range(count):
				x0[i] += dlon[i] * dt
				y0[i] += dlat[i] * dt
				z0[i] += dz[i] * dt
			
			# Transform latlon to target proj
			x, y, z = pyproj.transform(self.projGeog, self.projTo, x0, y0, z0, True)

		else:

			# Reproject to the dest coordinates
			x, y, z = pyproj.transform(self.projECEF, self.projTo, x, y, z, True)

		bounds = [FLT_MAX, FLT_MIN, FLT_MAX, FLT_MIN, FLT_MAX, FLT_MIN]
		# Expand the bounds for the new header.
		for i in range(count):
			if x[i] < bounds[0]: bounds[0] = x[i]
			if x[i] > bounds[1]: bounds[1] = x[i]
			if y[i] < bounds[2]: bounds[2] = y[i]
			if y[i] > bounds[3]: bounds[3] = y[i]
			if z[i] < bounds[4]: bounds[4] = z[i]
			if z[i] > bounds[5]: bounds[5] = z[i]

		return (x, y, z, bounds)

	def epochTransform(self, x, y, z, dt):
		'''
		Transform the coordinate using the procedure listed in Craymer (2006).
		x, y, z -- 	The coordinate arrays.
		dt      --  The time delta in decimal years.
		'''
		Txt = self.params.tx + self.params.dtx * dt
		Tyt = self.params.ty + self.params.dty * dt
		Tzt = self.params.tz + self.params.dtz * dt
		Dst = self.params.d + self.params.dd * dt
		Rxt = _sec2rad(self.params.rx + self.params.drx * dt)
		Ryt = _sec2rad(self.params.ry + self.params.dry * dt)
		Rzt = _sec2rad(self.params.rz + self.params.drz * dt)

		Dst += 1.

		for i in range(len(x)):
			x[i] = Txt + (Dst * x[i]) + (-Rzt * y[i]) + (Ryt * z[i])
			y[i] = Tyt + (Rzt * x[i]) + (Dst * y[i]) + (-Rxt * z[i])
			z[i] = Tzt + (-Ryt * x[i]) + (Rxt * y[i]) + 	(Dst * z[i])

	def getEllipsProps(self, proj):
		'''
		Return the ellipsoid parameters a and e2.
		'''
		if not self.a:
			# ellipsoid, axis, flattening, ref
			p = subprocess.Popen(['proj', '-le'], stdout=subprocess.PIPE)
			stdout, stderr = p.communicate()
			ellps = {}
			for x in [x for x in stdout.split('\n')]:
				x = tuple(x.split())
				try:
					ellps[x[0]] = x
				except: pass
			f = 1. / float(ellps['GRS80'][2].split('=')[1])
			self.a = float(ellps['GRS80'][1].split('=')[1])
			self.e2 = 2 * f - f * f
		return (self.a, self.e2)

	def initProjections(self):
		'''
		Initialize projections.
		'''
		if not self.params.fromSRS or not self.params.toSRS:
			raise Exception("SRSes are not set.")
		self.projFrom = pyproj.Proj(self.params.fromSRS)
		self.projTo = pyproj.Proj(self.params.toSRS)
		self.projECEF = pyproj.Proj("+proj=geocent +ellps=GRS80 +units=m +no_defs")
		self.projGeog = pyproj.Proj("+proj=latlon +ellps=GRS80 +no_defs")

	def loadHelmert(self):
		'''
		Load the transformation database.
		'''
		found = False
		with open(ITRF_FILE, "rU") as f:
			line = f.readline()
			while line:
				line = line.strip().split()
				if line[0] == self.ffrom:
					ffrom, fto = line[:2]
					epoch, tx, ty, tz, rx, ry, rz, d, dtx, dty, dtz, drx, dry, drz, dd = map(float, line[2:])
					self.params.epoch = epoch
					self.params.d = d / 1000000000. 		# Listed as ppb.
					self.params.dd = dd / 1000000000.
					self.params.tx = tx
					self.params.ty = ty
					self.params.tz = tz
					self.params.rx = rx
					self.params.ry = ry
					self.params.rz = rz
					self.params.dtx = dtx
					self.params.dty = dty 
					self.params.dtz = dtz
					self.params.drx = drx
					self.params.dry = dry
					self.params.drz = drz
					found = True
					break
				line = f.readline()

		if not found:
			raise Exception("Failed to find a transformation matching the parameters.")

if __name__ == '__main__':

	try:

		params = sys.argv[1:]
		ffrom = params[0]
		efrom = float(params[1])
		eto = float(params[2])
		fromSRS = params[3]
		toSRS = params[4]
		coords = map(float, params[5:])

		t = Transformer(ffrom, efrom, eto, fromSRS, toSRS)
		x = []
		y = []
		z = []
		for i in range(0, len(coords), 3):
			x.append(coords[i])
			y.append(coords[i+1])
			z.append(coords[i+2])

		x, y, z, bounds = t.transformPoints(x, y, z)
		#print x, y, z, bounds

	except:
		import traceback
		traceback.print_exc()
		print "Usage: nad83csrs.py <from frame> <from epoch> <to epoch> <from SRS> <to SRS> <coords [coords ...]>"
		print "  This program transforms coordinates from some reference frame/epoch to NAD83(CSRS) at another epoch."
		print "  from frame    -- The source's reference frame."
		print "  from epoch    -- The source epoch as a decimal."
		print "  to epoch      -- The destination epoch."
		print "  from SRS      -- The source's SRS; a proj string. Can use +init=epsg:0000 format for EPSG codes."
		print "  to SRS        -- The destination SRS."
		print "  coords        -- A coordinate in the form 'x y z'; can repeat as many times as allowed by the OS."
