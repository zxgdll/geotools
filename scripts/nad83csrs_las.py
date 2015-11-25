#!/usr/bin/env python
# -*- coding: utf-8 -*-

import laspy
import nad83csrs
import sys
import os

def transform_las(srcfile, dstdir, ffrom, efrom, eto, zone):
	'''
	Transforms a LAS file, or directory of LAS files, from one reference frame at 
	one epoch, to NAD83(CSRS) at another epoch.

	srcfile 	-- The source file or directory containing LAS files.
	dstdir		-- An output directory. Will be created if necessary.
	ffrom		-- The origin reference frame. Will be something like 'itrf90'.
	efrom 		-- The origin epoch. A decimal year like 1998.3.
	eto 		-- The destination epoch.
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

	i = 1
	for f in files:
		print 'Processing %u of %u: %s' % (i, len(files), f)
		src = laspy.file.File(f, mode = 'r')
		x = src.x
		y = src.y
		z = src.z
		nad83csrs.transform(x, y, z, ffrom, efrom, eto, zone)
		dst = laspy.file.File(os.path.join(dstdir, os.path.basename(f)), 
			header = src.header, vlrs = src.header.vlrs, mode = 'w')
		dst.points = src.points
		dst.x = x
		dst.y = y
		dst.z = z
		dst.close()
		src.close()
		i = i + 1

if __name__ == '__main__':

	try:
		srcfile, dstdir, ffrom, efrom, eto = sys.argv[1:6]
		zone = 0
		try:
			zone = int(sys.argv[6])
		except: pass
		nad83csrs.init()
		transform_las(srcfile, dstdir, ffrom, float(efrom), float(eto), zone)
		nad83csrs.cleanup()
	except Exception, e:
		import traceback
		traceback.print_exc()
		print 'Usage: nad83csrs_las.py <srcfile (or dir)> <dstdir> <origin reference frame> <origin epoch> <destination epoch> <zone>'
