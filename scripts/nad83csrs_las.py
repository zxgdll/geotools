#!/usr/bin/env python

import laspy
import nad83csrs
import sys
import os

nad83csrs.SHIFT_FILE = '/home/rob/Documents/gvb/NAD83v6VG.tif'

def transform_las(srcfile, dstdir, ffrom, efrom, eto, type, zone):
	
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
		for f in [x for x in os.listdir(srcfile) if x.lower().endswith('.las')]:
			files.append(os.path.join(srcfile, f))
	if not len(files):
		raise Exception('No files found in ' + str(srcfile))

	for f in files:

		src = laspy.file.File(f, mode = 'r')
		x = src.x
		y = src.y
		z = src.z
		print 'before', x
		x, y, z = nad83csrs.transform(x, y, z, ffrom, efrom, eto, type, zone)
		print 'after', x
		dst = laspy.file.File(os.path.join(dstdir, os.path.basename(f)), 
			header = src.header, vlrs = src.header.vlrs, mode = 'w')
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