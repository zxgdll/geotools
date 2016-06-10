#!/usr/bin/env python

import gdal
import os
import sys

def fmt(v):
	return '{:.6f}'.format(v)

ds = gdal.Open(sys.argv[1])
band = ds.GetRasterBand(1)
for r in range(ds.RasterYSize):
	print ' '.join(map(fmt, band.ReadAsArray(0, r, ds.RasterXSize, 1)[0,]))
