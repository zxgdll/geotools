#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script attempts to find an adjustment plane to fit two 3-dimensional
# datasets together.

import sys
import os
import itertools
import numpy as np
from osgeo import gdal
from scipy.optimize import curve_fit
from random import random

def grid_open(filename, band = 1):
	'''
	Open a raster and return the data array, the transform object
	the dataset object, the nodata value, the width and the height.
	'''
	ds = gdal.Open(filename)
	band = ds.GetRasterBand(1)
	trans = ds.GetGeoTransform()
	w = ds.RasterXSize * trans[1]
	h = ds.RasterYSize * -trans[5]
	return (band.ReadAsArray(), trans, ds, band.GetNoDataValue(), w, h)

def random_point(minx, miny, maxx, maxy):
	'''
	Return a random point within the given bounds.
	'''
	x = minx + random() * (maxx - minx)
	y = miny + random() * (maxy - miny)
	return x, y

def create_grid(popt, outfile, trans, w, h):
	'''
	Create a raster with the given transform and dimensions,
	whose elevations are determined by the plane parameters
	contained in the first argument.
	'''
	cols = int(w / trans[1]) + 1
	rows = int(h / -trans[5]) + 1
	print cols, rows
	drv = gdal.GetDriverByName('GTiff')
	ds = drv.Create(outfile, cols, rows, 1, gdal.GDT_Float32)
	ds.SetGeoTransform(trans)

	band = np.array([0.] * cols * rows, dtype='float32')
	band.shape = (rows, cols)

	for r in range(rows):
		for c in range(cols):
			x = c * trans[1] + trans[0]
			y = r * trans[5] + trans[3]
			z = x * popt[0] + y * popt[1] + popt[2]
			#print x, y, z, c, r
			band.itemset((r, c), z)
	ds.GetRasterBand(1).WriteArray(band)
	ds = None

def create_grid2(xdata, ydata, zdata, outfile, trans, w, h):
	'''
	Create a raster with the given transform and dimensions,
	whose elevations are determined by the plane parameters
	contained in the first argument.
	'''
	cols = int(w / trans[1]) + 1
	rows = int(h / -trans[5]) + 1
	print cols, rows
	drv = gdal.GetDriverByName('GTiff')
	ds = drv.Create(outfile, cols, rows, 1, gdal.GDT_Float32)
	ds.SetGeoTransform(trans)

	band = np.array([0.] * cols * rows, dtype='float32')
	band.shape = (rows, cols)

	for r in range(rows-1):
		for c in range(cols-1):
			x = xdata[rows - 2 - r,c]# * cols + c] #c * trans[1] + trans[0]
			y = ydata[rows - 2 - r,c]# * cols + c] #r * trans[5] + trans[3]
			z = zdata[rows - 2 - r,c]# * cols + c] #x * popt[0] + y * popt[1] + popt[2]
			#print x, y, z, c, r
			band.itemset((r, c), z)
	ds.GetRasterBand(1).WriteArray(band)
	ds = None

def grid_fit(file1, file2, mask, outfile, samples=100):
	'''
	Find the best-fit parameters for fitting the second filel
	to the first, and create an adjustment raster called outfile.
	The number of samples dictates the number of samples used to 
	find residuals.
	'''
	band1, trans1, ds1, nd1, w1, h1 = grid_open(file1)
	band2, trans2, ds2, nd2, w2, h2 = grid_open(file2)
	try:
		bandm, transm, dsm, ndm, wm, hm = grid_open(mask)
	except:
		print 'No mask.'

	# Use the minimum bounds for collecting samples.
	minx = max(trans1[0], trans2[0])
	maxx = min(trans1[0] + ds1.RasterXSize * trans1[1], trans2[0] + ds2.RasterXSize * trans2[1])
	maxy = min(trans1[3], trans2[3])
	miny = max(trans1[3] + ds1.RasterYSize * trans1[5], trans2[3] + ds2.RasterYSize * trans2[5])

	xydata = []
	zdata = []
	points = None
	if isinstance(samples, (float, int)):
		with open('fit_pts.csv', 'w') as wf:
			while len(zdata) < samples:
				pt = random_point(minx, miny, maxx, maxy)
				c1 = int(round((pt[0] - trans1[0]) / w1 * band1.shape[1]))
				r1 = int(round((trans1[3] - pt[1]) / h1 * band1.shape[0]))
				c2 = int(round((pt[0] - trans2[0]) / w2 * band2.shape[1]))
				r2 = int(round((trans2[3] - pt[1]) / h2 * band2.shape[0]))
				cm = int(round((pt[0] - transm[0]) / wm * bandm.shape[1]))
				rm = int(round((transm[3] - pt[1]) / hm * bandm.shape[0]))
				try:
					if bandm.item(rm, cm) > 0:
						z1 = band1.item(r1, c1)
						z2 = band2.item(r2, c2)
						if z1 != nd1 and z2 != nd2:
							xydata.append(pt)
							zdata.append(z1 - z2)
							wf.write('%f,%f,%f,%f\n' % (pt[0], pt[1], z1, z2))
				except Exception, e:
					print e
	else:
		raise Exception('Not implemented')


	xydata = np.array(xydata)
	zdata = np.array(zdata)

	def _f0(pt, a, b, c):
		return pt[:,0] * a + pt[:,1] * b + c

	def _f2(pt, a, b, c, d, e, f, g):
		return a * pt[:,0] ** 2 + b * pt[:,1] ** 2 + c * p[:,0] * p[:,1] + d * p[:,0] + e * p[:,1] + f + g

	popt, pcov = curve_fit(_f2, xydata, zdata)

	# Use the max bounds for building the gradient.
	minx = min(trans1[0], trans2[0])
	maxx = max(trans1[0] + ds1.RasterXSize * trans1[1], trans2[0] + ds2.RasterXSize * trans2[1])
	maxy = max(trans1[3], trans2[3])
	miny = min(trans1[3] + ds1.RasterYSize * trans1[5], trans2[3] + ds2.RasterYSize * trans2[5])

	create_grid(popt, outfile, [minx, 100., 0., maxy, 0., -100.], maxx-minx, maxy-miny)

	print popt
	print pcov

def poly_matrix(x, y, order=2):
    """ generate Matrix use with lstsq """
    import itertools
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i, j) in enumerate(ij):
        G[:, k] = x**i * y**j
    return G

def point_fit(file, outfile, samples):
	'''
	Find the best-fit parameters for fitting a plane to the point-set
	Create an adjustment raster called outfile.
	File is used to create the bounds and properties of the output.
	'''
	band1, trans1, ds1, nd1, w1, h1 = grid_open(file)

	# Use the minimum bounds for collecting samples.
	minx = trans1[0]
	maxx = trans1[0] + ds1.RasterXSize * trans1[1]
	maxy = trans1[3]
	miny = trans1[3] + ds1.RasterYSize * trans1[5]

	xydata = []
	xdata = []
	ydata = []
	zdata = []
	points = None

	for x, y, z in samples:
		xydata.append((x, y))
		xdata.append(x)
		ydata.append(y)
		zdata.append(z)

	xydata = np.array(xydata)
	xdata = np.array(xdata)
	ydata = np.array(ydata)
	zdata = np.array(zdata)

	def _f(pt, a, b, c):
		return pt[:,0] * a + pt[:,1] * b + c

	def _f2(pt, a, b, c, d, e, f, g):
		return a * pt[:,0] ** 2 + b * pt[:,1] ** 2 + c * pt[:,0] * pt[:,1] + d * pt[:,0] + e * pt[:,1] + f + g
	
	def polyfit2d(x, y, z, order=3):
	    ncols = (order + 1)**2
	    G = np.zeros((x.size, ncols))
	    ij = itertools.product(range(order+1), range(order+1))
	    for k, (i,j) in enumerate(ij):
	        G[:,k] = x**i * y**j
	    m, _, _, _ = np.linalg.lstsq(G, z)
	    return m

	def polyval2d(x, y, m):
	    order = int(np.sqrt(len(m))) - 1
	    ij = itertools.product(range(order+1), range(order+1))
	    z = np.zeros_like(x)
	    for a, (i,j) in zip(m, ij):
	        z += a * x**i * y**j
	    return z

	m = polyfit2d(xdata, ydata, zdata)

	# Evaluate it on a grid...
	nx, ny = ds1.RasterXSize, ds1.RasterYSize
	xx, yy = np.meshgrid(np.linspace(xdata.min(), xdata.max(), nx), np.linspace(ydata.min(), ydata.max(), ny))
	zz = polyval2d(xx, yy, m)

	create_grid2(xx, yy, zz, outfile, [minx, trans1[1], 0., maxy, 0., trans1[5]], maxx-minx, maxy-miny)


def load_pg_points(dbname, user, password, query):
	import psycopg2
	conn = psycopg2.connect("dbname=%s user=%s" % (dbname, user))
	cur = conn.cursor()
	points = []
	cur.execute(query)
	for x, y, z in cur:
		points.append((x, y, z))
	cur.close()
	conn.close()
	return points

if __name__ == '__main__':

	#try:
	import getopt

	opts, args = getopt.getopt(sys.argv[1:], 'f:g:m:o:p:q:u:s:d:')
	dbname = None
	dbuser = None
	dbpass = None
	query = None
	file1 = None
	file2 = None
	outfile = None
	samples = None
	mask = None
	for opt, val in opts:
		if opt == '-f':
			file1 = val
		elif opt == '-g':
			file2 = val
		elif opt == '-m':
			mask = val
		elif opt == '-u':
			dbuser = val
		elif opt == '-d':
			dbname = val
		elif opt == '-q':
			query = val
		elif opt == '-o':
			outfile = val
		elif opt == '-s':
			samples = int(samples)

	if dbname and query and outfile and file1:
		samples = load_pg_points(dbname, dbuser, dbpass, query)
		point_fit(file1, outfile, samples)
	elif file1 and file2 and outfile and mask and samples:
		grid_fit(file1, file2, mask, outfile, samples)
	else:
		raise Exception('Fail.')
	#except Exception, e:
	#	print e
	#	print 'Usage: grid_fit.py <file1>'
