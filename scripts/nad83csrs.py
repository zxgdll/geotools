#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
This script converts coordinates from and projection and any ITRF reference frame
to NAD83(CSRS), including the grid shifts to account for tectonic rotation and isostatic 
rebound. The script uses the grid shift raster published by NRCAN and transformation 
tables published for the ITRS.

The important method is to_nad83csrs, which takes these arguments:
	
	x, y, z 	-- 	3D Coordinate, scalar or array input is allowed.
	ffrom 		-- 	The source reference frame. This will be something like 'itrf90'.
	efrom		-- 	The source epoch. This is in decimal years.
	eto 		-- 	The destination epoch. The coordinates will be transformed to 
					NAD83(CSRS) at the 1997.0 epoch, then shifted to the desired epoch.

Notes from the IERS site:

	Note : These parameters are derived from those already published in the IERS
	Technical Notes indicated in the table above. The transformation parameters 
	should be used with the standard model (1) given below and are valid at the
	indicated epoch. 

	: XS :    : X :   : T1 :   :  D   -R3   R2 : : X :
	:    :    :   :   :    :   :               : :   :
	: YS :  = : Y : + : T2 : + :  R3   D   -R1 : : Y :                       (1)
	:    :    :   :   :    :   :               : :   :
	: ZS :    : Z :   : T3 :   : -R2   R1   D  : : Z :

	Where X,Y,Z are the coordinates in ITRF2000 and XS,YS,ZS are the coordinates in 
	the other frames.
	    
	On the other hand, for a given parameter P, its value at any epoch t 
	is obtained by using equation (2).
	                  
	P(t) = P(EPOCH) + P * (t - EPOCH)                                        (2)
	                                                          
	where EPOCH is the epoch indicated in the above table and P is the rate
	of that parameter.
'''

import sys
import numpy as np
from pyproj import Proj, transform as project
from scipy.interpolate import RectBivariateSpline
from osgeo import gdal

'''
This is the default path to the shift file. It can be overridden.
'''
SHIFT_FILE = 'NAD83v6VG.tif'

'''
This is the default path to the ITRF transform parameter file.
'''
ITRF_FILE = 'itrf.csv'

'''
This is the local datastore. In a database script, this will be replaced with 
a reference to DS.
'''
DATA = {}

def load_helmert():
	'''
	Loads the Helmert parameters from the local file into the data store.
	'''
	helmert_params = {}
	with open(ITRF_FILE, 'rU') as f:
		f.readline()
		line = f.readline()
		while line:
			try:
				cols = map(str.strip, line.split(','))
				frm, to = map(str.lower, cols[:2])
				epoch, tx, ty, tz, rx, ry, rz, ds, dtx, dty, dtz, drx, dry, drz, dds = map(float, cols[2:])
				helmert_params[frm] = (epoch, (tx, ty, tz, rx, ry, rz, ds), (dtx, dty, dtz, drx, dry, drz, dds))
			except Exception, e:
				pass
			line = f.readline()
	DATA['helmert_params'] = helmert_params


def load_shifts():
	'''
	Loads the shifts from the 3d shift grid, and initializes a set
	of interpolation functions for each dimension. These are stored
	in the data store.
	'''	
	ds = gdal.Open(SHIFT_FILE)
	trans = ds.GetGeoTransform()
	xcoords = np.arange(0, ds.RasterXSize)
	ycoords = np.arange(0, ds.RasterYSize)
	xs = ds.GetRasterBand(1).ReadAsArray().transpose()
	ys = ds.GetRasterBand(2).ReadAsArray().transpose()
	zs = ds.GetRasterBand(3).ReadAsArray().transpose()
	xspline = RectBivariateSpline(xcoords, ycoords, xs, kx = 1, ky = 1)
	yspline = RectBivariateSpline(xcoords, ycoords, ys, kx = 1, ky = 1)
	zspline = RectBivariateSpline(xcoords, ycoords, zs, kx = 1, ky = 1)
	DATA['shift_grids'] = (xspline, yspline, zspline)
	DATA['shift_trans'] = trans


def shift(x, y, efrom, eto):
	'''
	Returns the delta in each dimension (in mm) for the given range of 
	epochs.

	x 		-- The x coordinate(s).
	y 		-- The y coordinate(s).
	efrom 	-- The start epoch.
	eto 	-- The end epoch.
	'''
	sx, sy, sz = DATA['shift_grids']
	trans = DATA['shift_trans']
	dt = (eto - efrom)
	dx = [0.] * len(x)
	dy = [0.] * len(x)
	dz = [0.] * len(x)
	def _x(x, y, i):
		x = (x - trans[0]) / trans[1]
		y = (y - trans[3]) / trans[5]
		dx[i] = sx(x, y) / 1000. * dt
		dy[i] = sy(x, y) / 1000. * dt
		dz[i] = sz(x, y) / 1000. * dt
	map(_x, x, y, range(len(x)))
	return (dx, dy, dz)


def init():
	'''
	Initializes the data required for transformations.
	If cleanup is not called after, this creates a memory leak.
	'''
	load_helmert()
	load_shifts()


def cleanup():
	'''
	Free the stored data for GC.
	'''
	DATA['helmert_params'] = None
	DATA['shift_grids'] = None
	DATA['shift_trans'] = None


def epoch_transform(x, y, z, params, epoch):
	'''
	Transform the coordinate using the procedure listed in Craymer (2006).

	x 		-- 	The x coordinate(s).
	y 		-- 	The y coordinate(s).
	z 		-- 	The z coordinate(s).
	params 	-- 	The transformation parameters -- a tuple containing two tuples. One 
				containing the shift parameters (tx, ty, tz, ds) and one containing 
				the rates (dtx, dty, dtz, dds).
	epoch 	-- 	The epoch of the original coordinates.
	'''
	dt = epoch - params[0]
	tx, ty, tz, rx, ry, rz, d = params[1]
	dtx, dty, dtz, drx, dry, drz, dd = params[2]
	d = d / 1000000000. # Listed as ppb.
	dd = dd / 1000000000.
	a0 = tx + dtx * dt
	a1 = ty + dty * dt
	a2 = tz + dtz * dt
	bsx = 1. + d + dd * dt
	b01 = -rad(rz + drz * dt)
	b02 = rad(ry + dry * dt)
	b10 = rad(rz + drz * dt)
	b12 = -rad(rx + drx * dt)
	b20 = -rad(ry + dry * dt)
	b21 = rad(rx + drx * dt)
	def _trans(i):
		x[i] = a0 + bsx * x[i] + b01 * y[i] + b02 * z[i]
		y[i] = a1 + b10 * x[i] + bsx * y[i] + b12 * z[i]
		z[i] = a2 + b20 * x[i] + b21 * y[i] + bsx * z[i]
	map(_trans, range(len(x)))


def rad(arcsec):
	'''
	Convert the argument from arcseconds to radians.
	'''
	return arcsec * 4.84813681 / 1000000000.


def get_utm_srid(zone):
	'''
	Return the SRID for the WGS84 version of UTM, given the zone.
	'''
	return 32600 + int(zone)


def get_csrs_srid(zone, y = None):
	'''
	Figure out and return the SRID for a given NAD83(CSRS) UTM 
	projection given the zone, or the coordinate.

	zone	-- The x coordinate or zone, if known.
	y 		-- The y coordinate. If the first argument is zone, this must be None.
	'''
	if not y == None:
		x = float(zone)
		y = float(y)
		zone = int((x + 180.) / 6.) + 1;
		return get_csrs_srid(zone)
	else:
		zone = int(zone)
		if zone >= 11 and zone <= 13:
			return 2944 + zone
		elif zone >= 7 and zone <=10:
			return 3147 + zone
		elif zone >= 14 and zone <= 16:
			return 3144 + zone
		elif zone >= 17 and zone <= 21:
			return 2941 + zone
		else:
			raise Exception('Invalid zone for NAD83(CSRS): %u.' % (zone,))

def find_match(x0, y0, z0, x1, y1, z1, eto, zone = 0):
	'''
	x0, y0, z0 -- The points we know, in CSRS.
	x1, y1, z1 -- The input points. These are the ones we don't know the RF for.
	eto -- The epoch of the known points.
	zone -- The UTM zone or 0.
	'''
	results = None
	d = 99999.
	for e in range(19880, 20150):
		e = e  /10.
		for ffrom in DATA['helmert_params'].keys():
			if ffrom == 'nad83csrs': continue
			x, y, z = transform(x1, y1, z1, ffrom, e, eto, zone)
			dx0 = x0 - x
			dy0 = y0 - y
			dz0 = z0 - z
			d0 = abs(dx0) + abs(dy0) + abs(dz0)
			if d0 < d:
				result = (ffrom, e, dx0, dy0, dz0, x0, y0, z0, x1, y1, z1, x, y, z)
				d = d0
	return result

def transform(x, y, z, ffrom, efrom, eto, zone = 0):
	'''
	Transforms coordinate(s) from one available reference frame to NAD83(CSRS).
	x 		-- The x coordinate list.
	y 		-- The y coordinate list.
	z 		-- The z coordinate list.
	ffrom 	-- The starting reference frame (e.g. 'ITRF92')
	efrom 	-- The starting epoch in decimal years (e.g. 1995.2)
	eto 	-- The ending epoch in decimal years (e.g. 2013.8)
	zone	-- The origin UTM zone if the CRS is UTM. Otherwise latlon is assumed.
	'''	
	# If scalars are given, wrap them.
	useScalar = False
	if isinstance(x, str):
		x = [float(x)]
		y = [float(y)]
		z = [float(z)]
		useScalar = True
	elif not hasattr(x, '__len__'):
		x = [x]
		y = [y]
		z = [z]
		useScalar = True

	# Get the transformation params for the ref frames/epochs.
	params = DATA['helmert_params'][ffrom]

	# Get srids for to and from CRSes.
	if zone <= 0:
		type = 'latlon'
		fsrid = 4326
		tsrid = get_csrs_srid(x[0], y[0])
	else:
		type = 'nad83'
	 	fsrid = get_utm_srid(zone)
		tsrid = get_csrs_srid(zone)

	# Initialize projections.
	p1 = Proj('+init=EPSG:%u' % (fsrid,))
	p2 = Proj('+init=EPSG:4978')
	p3 = Proj('+init=EPSG:%u' % (tsrid,))
	p4 = Proj('+init=EPSG:4326')

	# Project to Cartesian 3D.
	# TODO: Unfortunately, this requires a copy
	x0, y0, z0 = project(p1, p2, x, y, z)

	# Transform using Helmert (etc.) params.	
	epoch_transform(x0, y0, z0, params, efrom)

	# Only use the grid shift if the epoch changes.
	dx, dy, dz = ([0.], [0.], [0.])
	if efrom != eto:
		x1, y1, z1 = project(p2, p4, x0, y0, z0)
		dx, dy, dz = shift(x1, y1, efrom, eto)

	# Reproject to the dest coordinates and add the shifts.
	x0, y0, z0 = project(p2, p3, x0, y0, z0)

	# If the epoch is changing, update with the shift values.
	# Otherwise, just write back to the input arrays.
	if efrom != eto:
		for i in range(len(x)):
			x[i] = x0[i] + dx[i]
			y[i] = y0[i] + dy[i]
			z[i] = z0[i] + dz[i]
	else:
		for i in range(len(x)):
			x[i] = x0[i]
			y[i] = y0[i]
			z[i] = z0[i]

	# Return values (not necessary if arrays passed in.)
	if useScalar:
		return x[0], y[0], z[0]
	else:
		return x, y, z

if __name__ == '__main__':

	try:
		argv = sys.argv[1:]
		mode = argv[0]
		if mode == 'transform':
			x, y, z = map(float, argv[1:4])
			ffrom = argv[4]
			efrom, eto = map(float, argv[5:7])
			zone = 0
			try:
				zone = int(argv[7])
			except: pass
			init()
			print transform(x, y, z, ffrom, efrom, eto, zone)
			cleanup()
		elif mode == 'find_match':
			x0, y0, z0 = map(float, argv[1:4])
			x1, y1, z1 = map(float, argv[4:7])
			eto = float(argv[7])
			zone = int(argv[8])
			init()
			print find_match(x0, y0, z0, x1, y1, z1, eto, zone)
			cleanup()
		else:
			raise Exception('Invalid mode:' + mode)
	except Exception, e:
		import traceback
		traceback.print_exc()
		print '''
Usage: nad83csrs.py <mode> [args]
	Available modes are: transform, find_match.
'''
