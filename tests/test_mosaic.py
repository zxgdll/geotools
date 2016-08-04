#!/usr/bin/env python

# Tests the mosaic program by performing the operation
# and checking the md5 sum of the output against
# a pre-determined value.

import os
import sys
import subprocess as sub

with open('data/mosaic.txt', 'r') as f:
	checksum1 = f.readline().strip()

sub.Popen(['../makefiles/mosaic', '-o', 'output/mosaic.tif', '-d', '10', 'data/mosaic_base.tif', 'data/mosaic_shapes.tif']).wait()
com = sub.Popen(['md5', '-q', 'output/mosaic.tif'], stderr=sub.PIPE, stdout=sub.PIPE)
checksum2 = com.stdout.readline().strip()

print 'Checksum match: ', checksum1 == checksum2
sys.exit(0 if checksum1 == checksum2 else 1)


