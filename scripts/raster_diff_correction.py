#!/usr/bin/env python

# This program computes the mean differences between valid pixels in a set of rasters.
# It builds a dependency graph by calculating the number of pixels of overlap between
# rasters and finding the maximum-cost path from each raster to a 'reference raster' --
# the one with no outgoing edges (i.e. the root). For any given raster, the adjustment 
# is the sum of the means of the differences along each edge from itself to the root.
# 
# This process works best if there is a node with no outgoing edges. 

import os
import sys
import gdal
import json
import numpy as np

def compute_stats(srcdir, nodata, cache = True):
	'''
	Computes the differences between rasters, and returns 
	a list containing tuples with:
		index1, index2, file1, file2, sum, count, mean
	
	If cache is set to false, the cache file (diffs.json) is loaded
	instead of calculating the actual differences. If the file is 
	not found, it is rebuilt.
	'''

	if cache:
		try:
			with open('diffs.json', 'rU') as f:
				return json.loads(f.read())
		except: pass

	files = sorted([x for x in os.listdir(srcdir) if x.endswith('.tif')])
	result = []

	for k in range(len(files)):
		for j in range(k + 1, len(files)):

			f = files[k]
			g = files[j]

			ds1 = gdal.Open(os.path.join(srcdir, f))
			ds2 = gdal.Open(os.path.join(srcdir, g))

			b1 = ds1.GetRasterBand(1)
			b2 = ds2.GetRasterBand(1)

			t1 = ds1.GetGeoTransform()
			t2 = ds2.GetGeoTransform()

			if t1[1] != t2[1] or t1[5] != t2[5]:
				raise Exception('Resolution mismatch.')

			xmin = max(t1[0], t2[0])
			xmax = min(t1[0] + ds1.RasterXSize * t1[1], t2[0] + ds2.RasterXSize * t2[1])
			ymin = max(t1[3] + ds1.RasterYSize * t1[5], t2[3] + ds2.RasterYSize * t2[5])
			ymax = min(t1[3], t2[3])

			c = 0
			v = 0.

			c1 = int((xmin - t1[0]) / t1[1])
			c2 = int((xmin - t2[0]) / t2[1])
			w = min(ds1.RasterXSize - c1, ds2.RasterXSize - c2)
			y = ymax

			# Read each row where rasters overlap
			while y > ymin:
				r1 = int((y - t1[3]) / t1[5])
				r2 = int((y - t2[3]) / t2[5])
				y += t1[5]

				v1 = b1.ReadAsArray(c1, r1, w, 1)
				v2 = b2.ReadAsArray(c2, r2, w, 1)

				# Build a mask to eliminate indices with at least one nodata.
				test = np.logical_and(np.logical_not(np.in1d(v1[0], [nodata])), np.logical_not(np.in1d(v2[0], [nodata])))

				# Difference the masked arrays.
				v0 = v1[0][test] - v2[0][test]

				if len(v0):
					v += v0.sum()
					c += len(v0)

			if c > 0:
				result.append((k, j, f, g, v, c, v / c))

	# Save the cache
	if cache:
		with open('diffs.json', 'w') as f:
			f.write(json.dumps(result))

	return result

class Edge(object):
	'''
	Represents an edge between two rasters (Nodes).
	'''

	def __init__(self, nfrom, nto, count, mean):
		self.nfrom = nfrom
		self.nto = nto
		self.count = count
		self.mean = mean

	def __repr__(self):
		return '[Edge: %s -> %s; count: %u; mean: %f]' % (self.nfrom, self.nto, self.count, self.mean)

	def __str__(self):
		return self.__repr__()


class Node(object):
	'''
	Represents a raster.
	'''
	def __init__(self, filename):
		self.filename = filename
		self.incoming = []
		self.outgoing = []

	def __repr__(self):
		return '[Node: %s; in %u, out %u]' % (self.filename, len(self.incoming), len(self.outgoing))

	def __str__(self):
		return self.__repr__()

def build_chains(pairs, reference):
	'''
	Builds and returns a series of 'chains' -- lists of
	edges from which shifts can be computed.
	'''
	
	# Build graph
	nodes = {}
	for pair in pairs:
		n1 = nodes.get(pair[2], None)
		if not n1:
			n1 = nodes[pair[2]] = Node(pair[2])
		n2 = nodes.get(pair[3], None)
		if not n2:
			n2 = nodes[pair[3]] = Node(pair[3])
		e1 = Edge(n1, n2, pair[5], -pair[6])
		e2 = Edge(n2, n1, pair[5], pair[6])
		n1.incoming.append(e1)
		n2.outgoing.append(e2)

	# Find the leaves
	leaves = []
	for k, n in nodes.iteritems():
		if not len(n.incoming):
			leaves.append(n)

	# Build the chains from the leaves
	chains = []
	for n in leaves:
		chain = []
		while n and len(n.outgoing):
			e = sorted(n.outgoing, key = lambda x: x.count, reverse = True)[0]
			n = e.nto
			chain.append(e)
		chains.append(chain)

	return chains

if __name__ == '__main__':

	try:
		srcdir = sys.argv[1]
		nodata = -9999.
		try:
			nodata = sys.argv[2]
		except: pass

		# Compute the differences
		result = compute_stats(srcdir, nodata)
		# Build the chains
		chains = build_chains(result, None)

		# Output a table of filenames and shifts.
		for c in chains:
			diff = 0.
			for e in c[::-1]:
				diff += e.mean
				print '%s,%f,%f' % (e.nfrom.filename, e.mean, diff)
	except Exception, e:
		print e
		print "Usage: diffs.py <src dir> [nodata value]"

