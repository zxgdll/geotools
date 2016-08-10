#!/usr/bin/env python

import sys
import os
import subprocess as sub

print 'Building las file.'
sub.Popen(['txt2las', '-a_srs', 'epsg:2956', '-i', 'data/lasgrid_data.txt', '-o', '/tmp/lasgrid_data.las']).wait()

print 'Loading checksums.'
args = {}
with open('data/lasgrid.txt', 'r') as f:
	line = f.readline()
	while line:
		try:
			line = map(str.strip, line.split())
			args[line[0]] = line[1]
		except: pass
		line = f.readline()

print 'Testing...'
for arg in args.keys():
	print 'Testing', arg
	sub.Popen(['../makefiles/lasgrid-app', '-v', '-o', '/tmp/lasgrid_test.tif', '-t', arg, '-r', '1', '-d', '0.77', '-s', '2956', '/tmp/lasgrid_data.las']).wait()
	com = sub.Popen(['md5sum', '/tmp/lasgrid_test.tif'], stdout=sub.PIPE, stderr=sub.PIPE)
	print com.stdout.read()
	if com.returncode == 0:
		if args[arg] != com.stdout.read().strip().split()[0].strip():
			print 'Failed on', arg
			sys.exit(1)
	else:
		print 'Failed'
		sys.exit(1)
print 'Done'
sys.exit(0)
