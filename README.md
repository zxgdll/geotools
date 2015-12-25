# geotools

Contains scripts and things that I've found useful.

## las2csrs

This program is experimental. Please do not use it for important work.

This is a c++ program that transforms a point cloud in LAS format from any ITRF reference frame/epoch to NAD83(CSRS) at any epoch. The algorithm is based on Craymer, M. R. (2006). The evolution of NAD83 in Canada, 60(2), 151–165.

## gvb2tif.py

This program converts NRCAN's NAD83(CSRS) velocity grid to a GeoTiff format. The grid contains horizontal and 
vertical velocities for calculating tectonic and isostatic change of coordinates over time.

![RGB Rendering of XYZ Bands](http://dijital.ca/img/gvb2.jpg "RGB Rendering of XYZ Bands")

