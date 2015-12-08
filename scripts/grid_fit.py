#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script attempts to find an adjustment plane to fit two 3-dimensional
# datasets together.

import sys
import os
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

	def _f(pt, a, b, c):
		return pt[:,0] * a + pt[:,1] * b + c

	popt, pcov = curve_fit(_f, xydata, zdata)

	# Use the max bounds for building the gradient.
	minx = min(trans1[0], trans2[0])
	maxx = max(trans1[0] + ds1.RasterXSize * trans1[1], trans2[0] + ds2.RasterXSize * trans2[1])
	maxy = max(trans1[3], trans2[3])
	miny = min(trans1[3] + ds1.RasterYSize * trans1[5], trans2[3] + ds2.RasterYSize * trans2[5])

	create_grid(popt, outfile, [minx, 100., 0., maxy, 0., -100.], maxx-minx, maxy-miny)

	print popt
	print pcov

def load_pg_points(dbname, user, query):
	import psycopg2
	conn = psycopg2.connect("dbname=%s user=%s" % (dbname, user))
	cur = conn.cursor()
	points = []
	for x, y in cur.execute(query):
		points.append((x, y))
	cur.close()
	conn.close()
	return points

if __name__ == '__main__':

	#try:
	file1 = sys.argv[1]
	file2 = sys.argv[2]
	mask = sys.argv[3]
	outfile = sys.argv[4]
	samples = sys.argv[5]
	if samples.startswith('pg'):
		dbname, user, query = sys.argv[6:9]
		samples = load_pg_points(dbname, user, query)
	else:
		samples = int(samples)
	grid_fit(file1, file2, mask, outfile, samples)
	#except Exception, e:
	#	print e
	#	print 'Usage: grid_fit.py <file1>'
