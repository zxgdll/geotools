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

# tx, ty, tz, scale, rx, ry, rz
# m, m, m, ppm, sec, sec, sec
# rates are /yr.
helmert_params = {}

def load_helmert():
	'''
	Loads the Helmert parameters from the local file.
	'''
	def remap_params():
		'''
		This method modifies the helmert param dict to provide reverse parameters
		between frames that aren't represented.
		'''
		def reverse(item):
			item = list(item)
			item[0] = tuple([-x for x in item[0]])
			item[1] = tuple([-x for x in item[1]])
			return tuple(item)

		for k0, v0 in list(helmert_params.iteritems()):
			# from, to
			for k1, v1 in v0.iteritems():
				# to, epoch
				d0 = helmert_params.get(k1, None)
				if not d0:
					d0 = helmert_params[k1] = {}
				d1 = helmert_params[k1].get(k0, None)
				if not d1:
					d1 = helmert_params[k1][k0] = {}
				for k2, v2 in v1.iteritems():
					#print k2, v2
					# epoch, params					
					if not d1.get(k2, None):
						d1[k2] = reverse(v2)

	with open('itrf.csv', 'rU') as f:
		f.readline()
		line = f.readline()
		while line:
			try:
				cols = map(str.strip, line.split(','))
				frm, to = map(str.lower, cols[:2])
				epoch, tx, ty, tz, d, rx, ry, rz, dtx, dty, dtz, dd, drx, dry, drz = map(float, cols[2:])
				if not helmert_params.get(frm, None):
					helmert_params[frm] = {}
				if not helmert_params[frm].get(to, None):
					helmert_params[frm][to] = {}
				helmert_params[frm][to][epoch] = ((tx, ty, tz, d, rx, ry, rz), (dtx, dty, dtz, dd, drx, dry, drz))
			except Exception, e:
				print 'error', e
				print line
			line = f.readline()

	remap_params()

def get_params(ffrom, efrom):
	'''
	Returns the transformation parameters for transformations between reference 
	frames at specific epochs.

	ffrom 	-- The starting reference frame (e.g. 'ITRF92')
	efrom 	-- The starting epoch in decimal years (e.g. 1995.2)
	'''
	# Search for the transformation path
	fto = 'itrf96'
	t = dict(helmert_params)
	n = t[ffrom]
	n['_parent'] = None
	n['_name'] = ffrom
	path = [n]
	found = False
	visited = set([ffrom])
	while path and not found:
		n = path.pop()
		for k, v in [x for x in t[n['_name']].iteritems() if not x[0].startswith('_') and not x[0] in visited]:
			visited.add(k)
			v = t[k]
			v['_parent'] = n
			v['_name'] = k
			path.append(v)
			if k == fto:
				found = True
				break

	n = helmert_params[k]
	path = []
	while n:
		path.append(n['_name'])
		n = n['_parent']
	path = path[::-1]

	epoch = efrom
	result = [0.] * 7
	for i in range(len(path)-1):
		a = path[i]
		b = path[i+1]
		t = helmert_params[a][b].keys()[0]
		params, rates = helmert_params[a][b][t]
		#print a, b, t, epoch, params, rates
		for j in range(7):
			result[j] = result[j] + params[j] + rates[j] * (epoch - t)
		epoch = t

	a = 'itrf96'
	b = 'nad83csrs'
	t = helmert_params[a][b].keys()[0]
	params, rates = helmert_params[a][b][t]
	#print a, b, t, params, rates
	for j in range(7):
		result[j] = result[j] + params[j]

	return tuple(result)

shift_funcs = None

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
	global shift_funcs
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
		shift_funcs = (fx, fy, fz)
	
	fx, fy, fz = shift_funcs
	dt = (eto - efrom)
	c = int((x - trans[0]) / trans[1])
	r = int((y - trans[3]) / trans[5])
	
	return (fx(x, y) / 1000. * dt, fy(x, y) / 1000. * dt, fz(x, y) / 1000. * dt)

def to_nad83csrs(x, y, z, sfrom, ffrom, efrom, eto):
	'''
	Transforms a coordinate from one reference frame to another using the proj library.
	x 		-- The x coordinate.
	y 		-- The y coordinate.
	z 		-- The z coordinate.
	sfrom 	-- The starting SRID.
	ffrom 	-- The starting reference frame (e.g. 'ITRF92')
	efrom 	-- The starting epoch in decimal years (e.g. 1995.2)
	eto 	-- The ending epoch in decimal years (e.g. 2013.8)
	'''	
	# Load the transform table
	load_helmert()

	# Get the transformation params for the ref frames/epochs.
	params = get_params(ffrom, efrom)
	# Reproject from the source coords to the target ref frame lat/lon.
	p1 = pyproj.Proj('+init=EPSG:%u +towgs84=%f,%f,%f,%f,%f,%f,%f' % tuple([sfrom] + list(params)))
	p2 = pyproj.Proj('+init=EPSG:4326')
	x0, y0, z0 = pyproj.transform(p1, p2, x, y, z)

	# Only use the grid shift if the epoch changes.
	if efrom != eto:
		dx, dy, dz = get_shifts(x0, y0, efrom, eto)
	else:
		dx, dy, dz = (0., 0., 0.)

	# Reproject to the dest coordinates and add the shifts.
	p3 = pyproj.Proj('+init=EPSG:2956')
	x1, y1, z1 = pyproj.transform(p1, p3, x, y, z)

	# Return the shifted coordinate.	
	return (x1 + dx, y1 + dy, z1 + dz)

print (470000., 6520000., 200.)
print to_nad83csrs(470000., 6520000., 200., 26912, 'itrf90', 1995., 2015.)

