#!/usr/bin/env python
# -*- coding: utf-8 -*-

import laspy
import nad83csrs
import sys
import os

def transform_las(srcfile, dstdir, ffrom, efrom, eto, type, zone):
	'''
	Transforms a LAS file, or directory of LAS files, from one reference frame at 
	one epoch, to NAD83(CSRS) at another epoch.

	srcfile 	-- The source file or directory containing LAS files.
	dstdir		-- An output directory. Will be created if necessary.
	ffrom		-- The origin reference frame. Will be something like 'itrf90'.
	efrom 		-- The origin epoch. A decimal year like 1998.3.
	eto 		-- The destination epoch.
	type		-- 'latlon' or 'nad83'.
	zone 		-- The UTM zone, if type is 'nad83'. Ignored otherwise.	
	'''

	if srcfile == dstdir:
		raise Exception('Destination and source are the same.')

	if os.path.isfile(dstdir):
		raise Exception(str(dstdir) + ' is a file.')

	if not os.path.exists(dstdir):
		try:
			os.makedirs(dstdir)
		except: pass
	
	files = []

	if os.path.isfile(srcfile):
		files.append(srcfile)
	elif os.path.isdir(srcfile):
		for f in [x for x in os.listdir(srcfile) if x.lower().endswith('.las') and not os.path.exists(os.path.join(dstdir, x))]:
			files.append(os.path.join(srcfile, f))
	
	if not len(files):
		raise Exception('No files found in ' + str(srcfile))

	for f in files:
		print 'Processing ', f
		src = laspy.file.File(f, mode = 'r')
		x = src.x
		y = src.y
		z = src.z
		nad83csrs.transform(x, y, z, ffrom, efrom, eto, type, zone)
		dst = laspy.file.File(os.path.join(dstdir, os.path.basename(f)), 
			header = src.header, vlrs = src.header.vlrs, mode = 'w')
		dst.points = src.points
		dst.x = x
		dst.y = y
		dst.z = z
		dst.close()
		src.close()

if __name__ == '__main__':

	try:
		srcfile, dstdir, ffrom, efrom, eto, type = sys.argv[1:7]
		if type == 'nad83':
			zone = sys.argv[7]
		nad83csrs.init()
		transform_las(srcfile, dstdir, ffrom, float(efrom), float(eto), type, zone)
		nad83csrs.cleanup()
	except Exception, e:
		import traceback
		traceback.print_exc()
		print 'Usage: nad83csrs_las.py <srcfile (or dir)> <dstdir> <origin reference frame> <origin epoch> <destination epoch> <type (latlon|nad83)> <zone>'
