#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
This script creates a 6-band GeoTiff from the 3-element velocity grid provided by NRCAN
for transforming spatial data in NAD83(CSRS) between various epochs. The file's CRS is WGS84
(EPSG:4326). The source grids are provided with the TRX program, available on the NRCAN
website and only for Windows. Once you have the grid files, you can process them on any 
system.

NOTE: This software is experimental and there are no guarantees, epxresed or implied. Use 
the official TRX program for all and any work that matters.

Download TRX here:

	http://www.nrcan.gc.ca/earth-sciences/geomatics/geodetic-reference-systems/tools-applications/10925

The band listing is:

	1, 2, 3: X, Y and Z shifts per annum in mm.
	4, 5, 6: Accuracy values for X, Y and Z.

For more information on the velocity grid, see: 
	
	Craymer, M. R., Henton, J. A., & Piraszewksi, M. (2009). Predicting Present-Day Rates of 
	Glacial Isostatic Adjustment Using a Smoothed GPS-Based Velocity Field for the 
	Reconciliation of NAD83 Reference Frames in Canada. In Workshop on Monitoring North American 
	Geoid Change (p. 1).

For more on the file format, see:

	Craymer, M. (2014). NTV2_3D Grid File Format.

Run this program like so:

	./gvb2tif.py NAD83v6VG.gvb NAD83v6VG.tif
'''

import struct
import os
import sys
import numpy as np
from osgeo import gdal, osr

try:
	src = sys.argv[1]
	dst = sys.argv[2]
except:
	print 'Usage: ./gvb2tif.py <gvb file> <tif file>'
	sys.exit(1)

header = {}

def read_name(f):
	'''
	Read a header name from the file.
	'''
	return ''.join(struct.unpack('c' * 8, f.read(8))).strip()

def read_head(f, h):
	'''
	Read a string header and set it on the header dict.
	'''
	name = read_name(f)
	value = struct.unpack('c' * 8, f.read(8))
	f.read(8)
	h[name] = ''.join(value).strip()

def read_ihead(f, h):
	'''
	Read an integer header and set it on the header dict.
	'''
	name = read_name(f)
	value = struct.unpack('i' * 4, f.read(16))
	h[name] = value[0]

def read_dhead(f, h):
	'''
	Read a double header and set it on the header dict.
	'''
	name = read_name(f)
	value = struct.unpack('d' * 2, f.read(16))
	h[name] = value[0]

with open(src, 'rb') as f:

	# Read the headers. Not really doing anything with the result.
	for i in range(3):
		read_ihead(f, header)
	for i in range(4):
		read_head(f, header)
	for i in range(4):
		read_dhead(f, header)
	for i in range(4):
		read_head(f, header)
	for i in range(6):
		read_dhead(f, header)
	read_ihead(f, header)

	# Compute dimensions.
	rows = int((header['N_LAT']-header['S_LAT']) / header['LAT_INC'] + 1)
	cols = int((header['W_LON']-header['E_LON']) / header['LON_INC'] + 1)

	# Load the data -- it's unformatted fortran, so this is eash. But it
	# needs to be reshaped and mirrored
	data = np.fliplr(np.flipud(np.reshape(np.fromfile(f, dtype = np.float32, count = rows * cols * 6), (rows, cols, 6))))

	# Create the TIF.
	drv = gdal.GetDriverByName('GTiff')
	ds = drv.Create(dst, cols, rows, 6, gdal.GDT_Float32)

	# Set Transform. The coordinates are given in arcsec.
	trans = [-header['W_LON'] / 3600.0, header['LON_INC'] / 3600.0, 0, header['N_LAT'] / 3600.0, 0, -header['LAT_INC'] / 3600.0]
	ds.SetGeoTransform(trans)

	# Set SRS
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(4326)
	ds.SetProjection(srs.ExportToWkt())

	# Write value bands.
	ds.GetRasterBand(1).WriteArray(data[...,0])
	ds.GetRasterBand(2).WriteArray(data[...,1])
	ds.GetRasterBand(3).WriteArray(data[...,2])
	# Write accuracy bands.
	ds.GetRasterBand(4).WriteArray(data[...,3])
	ds.GetRasterBand(5).WriteArray(data[...,4])
	ds.GetRasterBand(6).WriteArray(data[...,5])

print 'Done'



