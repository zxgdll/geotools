#!/usr/bin/env python

# Tests the mosaic program by performing the operation
# and checking the md5 sum of the output against
# a pre-determined value.

import os
import sys
import subprocess as sub

print 'Loading checksum.'
with open('data/mosaic.txt', 'r') as f:
	checksum1 = f.readline().strip()

print 'Building mosaic.'
cmd = ['../makefiles/mosaic-app', '-v', '-o', 'output/mosaic.tif', '-d', '10', 'data/mosaic_base.tif', 'data/mosaic_shapes.tif']
print ' '.join(cmd)
sub.Popen(cmd, env={"OMP_NUM_THREADS" : "2"}).wait()

print 'Computing checksum.'
com = sub.Popen(['md5sum', 'output/mosaic.tif'], stderr=sub.PIPE, stdout=sub.PIPE)
checksum2 = com.stdout.readline().strip().split()[0]

if checksum1 != checksum2:
	print 'Checksum mismatch'
	sys.exit(1)

print 'Done'
sys.exit(0)


