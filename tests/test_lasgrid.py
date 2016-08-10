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
	cmd = ['../makefiles/lasgrid-app', '-d', '-1', '-o', '/tmp/lasgrid_test.tif', '-t', arg, '-r', '1', '-s', '2956', '/tmp/lasgrid_data.las']
	print ' '.join(cmd)
	com = sub.Popen(cmd, stdout=sub.PIPE, stderr=sub.PIPE)
	com.wait()
	print com.stdout.read()
	print com.stderr.read()
	if com.returncode !=0:
		print 'Failed'
		sys.exit(1)
	com = sub.Popen(['md5sum', '/tmp/lasgrid_test.tif'], stdout=sub.PIPE, stderr=sub.PIPE)
	com.wait()
	stdout = com.stdout.read()
	stderr = com.stderr.read()
	print 'Result:', stdout
	if com.returncode == 0:
		md5 = stdout.strip().split()[0].strip()
		if args[arg] != md5:
			print 'Failed on', arg, md5, 'should be', args[arg]
			sys.exit(1)
	else:
		print 'Failed'
		sys.exit(1)
print 'Done'
sys.exit(0)
