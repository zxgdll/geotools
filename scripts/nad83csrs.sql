 CREATE EXTENSION IF NOT EXISTS plpythonu;


-- Transforms a PostGIS geometry in any coordinate reference system at any epoch and reference frame
-- to NAD83(CSRS) at any epoch.
-- 
-- geom 	A 3D point geometry in WKB form.
-- ffrom	The origin reference frame. This will be something like 'ITRF90'.
-- efrom	The origin epoch. Epochs are expressed as decimal years.
-- eto		The destination epoch.
-- fsrid    The SRID of the origin geometry.
-- tsrid    The SRID of the destination geometry.
CREATE OR REPLACE FUNCTION ToNAD83CSRS(geom bytea, ffrom text, efrom double precision, eto double precision, fsrid integer, tsrid integer)
RETURNS bytea
AS $$

import sys
import os

script_base = '/Users/robskelly/Documents/geotools'
script_path = script_base + '/scripts'
if not script_path in sys.path:
	sys.path.append(script_path)
	os.environ['PATH'] += os.pathsep + script_base + '/bin'

import nad83csrs
from shapely.geometry import Point
from shapely.wkb import loads, dumps

nad83csrs.SHIFT_FILE = script_base + '/share/NAD83v6VG.tif'
nad83csrs.ITRF_FILE = script_base + '/share/itrf.csv'

pt = loads(geom)
t = nad83csrs.Transformer(ffrom.lower(), efrom, eto, fsrid, tsrid)
x, y, z, bounds = t.transformPoints([pt.x], [pt.y], [pt.z])
return dumps(Point(x[0], y[0], z[0]))

$$ LANGUAGE 'plpythonu';


-- Transforms a PostGIS geometry in any coordinate reference system at any epoch and reference frame
-- to NAD83(CSRS) at any epoch.
-- 
-- geoms 	An array of 3D point geometries in WKB form.
-- ffrom	The origin reference frame. This will be something like 'ITRF90'.
-- efrom	The origin epoch. Epochs are expressed as decimal years.
-- eto		The destination epoch.
-- fsrid    The SRID of the origin geometry.
-- tsrid    The SRID of the destination geometry.
CREATE OR REPLACE FUNCTION ToNAD83CSRS(geoms bytea[], ffrom text, efrom double precision, eto double precision, fsrid integer, tsrid integer)
RETURNS bytea[]
AS $$

import sys

script_base = '/Users/robskelly/Documents/geotools'
script_path = script_base + '/scripts'
if not script_path in sys.path:
	sys.path.append(script_path)
	os.environ['PATH'] += os.pathsep + script_base + '/bin'

import nad83csrs
from shapely.geometry import Point
from shapely.wkb import loads, dumps

nad83csrs.SHIFT_FILE = script_base + '/share/NAD83v6VG.tif'
nad83csrs.ITRF_FILE = script_base + '/share/itrf.csv'

count = len(geoms)
x = [None] * count
y = [None] * count
z = [None] * count

for i in range(count):
	pt = loads(geoms[i])
	x[i] = pt.x
	y[i] = pt.y
	z[i] = pt.z

t = nad83csrs.Transformer(ffrom.lower(), efrom, eto, fsrid, tsrid)
x, y, z, bounds = t.transformPoints(x, y, z)

result = [None] * count
for i in range(count):
	result[i] = dumps(Point(x[i], y[i], z[i]))

return result

$$ LANGUAGE 'plpythonu';
