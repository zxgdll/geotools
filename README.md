# geotools

Contains scripts and things that I've found useful.

## las2csrs

This program is experimental. Please do not use it for important work.

This is a c++ program that transforms a point cloud in LAS format from any ITRF reference frame/epoch to NAD83(CSRS) at any epoch. The algorithm is based on Craymer, M. R. (2006). The evolution of NAD83 in Canada, 60(2), 151–165.

## lasgrid

This is a simple program that grids one or more LAS files into a GeoTiff. The user can select the classes, and a number of methods, such as mean, min, max, stddev, density, etc.

## lasstats

This program computes zonal statistics from a point cloud in LAS format from either a collection of polygons in Shapefile format, or from a classification raster. If a Shapefile is used, each polygon is updated with columns for the statistics. If a raster is used (it must be the Byte type), statistics are computed for each unique cell value and output as a CSV file.

The user may select which point classes are considered, but this behaviour is different depending on whether raster of vector is used. For rasters, statistics are computed for each cell value *and* for each point class. So, if there are two cell values and three point classes, there are six records in the table. For vectors, statistics are grouped for all given classes.

This program enables the user to develop something like this, where the x-axis is cell IDs, and the y-axis is canopy height.

[![Canopy height quartiles.](http://dijital.ca/files/geotools/lasstats.jpg](http://dijital.ca/files/geotools/lasstats.html)
Click on the image to view an interactive map.

## gvb2tif.py

This program converts NRCAN's NAD83(CSRS) velocity grid to a GeoTiff format. The grid contains horizontal and 
vertical velocities for calculating tectonic and isostatic change of coordinates over time.

![RGB Rendering of XYZ Bands](http://dijital.ca/img/gvb2.jpg "RGB Rendering of XYZ Bands")

