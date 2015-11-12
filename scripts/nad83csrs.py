#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
This script converts coordinates from and projection and any ITRF reference frame
to NAD83(CSRS), including the grid shifts to account for tectonic rotation and isostatic 
rebound. The script uses the grid shift raster published by NRCAN and transformation 
tables published for the ITRS.

The important method is to_nad83csrs, which takes these arguments:
	
	x, y, z 	-- 	3D Coordinate.
	sfrom	 	-- 	The EPSG ID of the original projection. This is also known as the SRID.
	ffrom 		-- 	The source reference frame. This will be something like 'itrf90'.
	efrom		-- 	The source epoch. This is in decimal years.
	eto 		-- 	The destination epoch. The coordinates will be transformed to 
					NAD83(CSRS) at the 1997.0 epoch, then shifted to the desired epoch.

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
        
                  .
P(t) = P(EPOCH) + P * (t - EPOCH)                                        (2)

                                                          .
where EPOCH is the epoch indicated in the above table and P is the rate
of that parameter.

'''

import pyproj
import sys
import numpy as np
from scipy.interpolate import interp2d
from pprint import pprint

try:
	DS
except:
	DS = {}

def load_helmert():
	'''
	Loads the Helmert parameters from the local file.
	'''
	helmert_params = DS.get('helmert_params')
	if not helmert_params:
		helmert_params = {}
		with open('itrf.csv', 'rU') as f:
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
		DS['params'] = helmert_params
	return helmert_params

def get_shifts(x, y, efrom, eto):
	'''
	Loads the shifts from the 3d shift grid, and initializes a set
	of interpolation functions for each dimension.

	Returns the delta in each dimension (in m) for the given range of 
	epochs.

	x 		-- The x coordinate.
	y 		-- The y coordinate.
	efrom 	-- The start epoch.
	eto 	-- The end epoch.
	'''
	shift_funcs = DS.get('shift_funcs', None)
	if not shift_funcs:
		# Load the raster and create interpolation functions
		# if they do not already exist.
		from osgeo import gdal
		ds = gdal.Open('NAD83v6VG.tif')
		trans = ds.GetGeoTransform()
		xx = np.arange(trans[0], trans[0] + ds.RasterXSize * trans[1], trans[1])
		yy = np.arange(trans[3], trans[3] + ds.RasterYSize * trans[5], trans[5])
		zz0 = ds.GetRasterBand(1).ReadAsArray()
		fx = interp2d(xx, yy, zz0, kind = 'cubic', bounds_error = True)
		zz1 = ds.GetRasterBand(2).ReadAsArray()
		fy = interp2d(xx, yy, zz1, kind = 'cubic', bounds_error = True)
		zz2 = ds.GetRasterBand(3).ReadAsArray()
		fz = interp2d(xx, yy, zz2, kind = 'cubic', bounds_error = True)
		shift_funcs = DS['shift_funcs'] = (fx, fy, fz)
	fx, fy, fz = shift_funcs
	dt = (eto - efrom)
	c = int((x - trans[0]) / trans[1])
	r = int((y - trans[3]) / trans[5])
	print 'shifts %0.9f %0.9f %0.9f' % (fx(x, y), fy(x, y), fz(x, y))
	return (fx(x, y) / 1000. * dt, fy(x, y) / 1000. * dt, fz(x, y) / 1000. * dt)

def transform(x, y, z, params, epoch):
	'''
	Transform the coordinate using the procedure listed in Craymer (2006).

	x 		-- 	The x coordinate.
	y 		-- 	The y coordinate.
	z 		-- 	The z coordinate.
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
	a = np.matrix([[tx + dtx * dt], [ty + dty * dt], [tz + dtz * dt]])
	#print 'params %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f' % (tx + dtx * dt, ty + dty * dt, tz + dtz * dt, d + dd * dt, rx + drx * dt, ry + dry * dt, rz + drz * dt)
	b = np.matrix([
		[1 + d + dd * dt, 		-rad(rz + drz * dt), 	rad(ry + dry * dt)],
		[rad(rz + drz * dt), 	1 + d + dd * dt, 		-rad(rx + drx * dt)],
		[-rad(ry + dry * dt), 	rad(rx + drx * dt),		1 + d + dd * dt]
	])
	c = np.matrix([[x], [y], [z]])
	d = np.add(a, np.dot(b, c))
	return (d[0,0], d[1,0], d[2,0])

def rad(arcsec):
	'''
	Convert the argument from arcseconds to radians.
	'''
	return arcsec * 4.84813681 / 1000000000.

def deg(arcsec):
	return arcsec * 0.0166667

def to_nad83csrs(x, y, z, ffrom, efrom, eto):
	'''
	Transforms a coordinate from one available reference frame to NAD83(CSRS).
	x 		-- The x coordinate.
	y 		-- The y coordinate.
	z 		-- The z coordinate.
	sfrom 	-- The starting SRID.
	ffrom 	-- The starting reference frame (e.g. 'ITRF92')
	efrom 	-- The starting epoch in decimal years (e.g. 1995.2)
	eto 	-- The ending epoch in decimal years (e.g. 2013.8)
	'''	
	# Load the transform table and get the transformation params for the ref frames/epochs.
	helmert_params = load_helmert()
	params = helmert_params[ffrom]

	# Reproject from the source coords to 3D Cartesian.
	p1 = pyproj.Proj('+init=EPSG:32612')
	p2 = pyproj.Proj('+init=EPSG:4978')
	x0, y0, z0 = pyproj.transform(p1, p2, x, y, z)
	#print '3D %0.9f %0.9f %0.9f' % (x0, y0, z0)

	# Transform using Helmert params.	
	x1, y1, z1 = transform(x0, y0, z0, params, efrom)
	#print 'trans %0.9f %0.9f %0.9f' % (x1, y1, z1,)
	# Only use the grid shift if the epoch changes.
	if efrom != eto:
		p4 = pyproj.Proj('+init=EPSG:4326')
		x3, y3, z3 = pyproj.transform(p2, p4, x1, y1, z1)
		dx, dy, dz = get_shifts(x3, y3, 1997., eto)
	else:
		dx, dy, dz = (0., 0., 0.)
	#print 'grid %0.9f %0.9f %0.9f' % (dx, dy, dz,)

	# Reproject to the dest coordinates and add the shifts.
	p3 = pyproj.Proj('+init=EPSG:32612')
	x2, y2, z2 = pyproj.transform(p2, p3, x1, y1, z1)
	#print 'final %0.9f %0.9f %0.9f' % (x2, y2, z2,)

	# Return the shifted coordinate.	
	return (x2 + dx, y2 + dy, z2 + dz)

print (470000., 6520000., 200.)
print to_nad83csrs(470000., 6520000., 200., 'itrf90', 1990., 2015.)

