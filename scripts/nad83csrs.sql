 CREATE EXTENSION IF NOT EXISTS plpythonu;


-- Transforms a PostGIS geometry in any coordinate reference system at any epoch and reference frame
-- to NAD83(CSRS) at any epoch.
-- 
-- geom 	A 3D point geometry in WKB form.
-- ffrom	The origin reference frame. This will be something like 'ITRF90'.
-- efrom	The origin epoch. Epochs are expressed as decimal years.
-- eto		The destination epoch.
-- fromsrs  The proj string of the origin geometry.
-- tosrs    The proj string of the destination geometry.
CREATE OR REPLACE FUNCTION _ToNAD83CSRS(geom bytea, ffrom text, efrom float8, eto float8, fromsrs text, tosrs text)
RETURNS bytea
AS $$

import sys
import os

script_base = '/Users/robskelly/Documents/geotools' # TODO: This needs to go.
script_path = script_base + '/scripts'
if not script_path in sys.path:
	sys.path.append(script_path)
	os.environ['PATH'] += os.pathsep + script_base + '/bin'
	os.environ['PROJ_LIB'] = os.pathsep + '/usr/local/share/proj'

import nad83csrs
from shapely.geometry import Point
from shapely.wkb import loads, dumps

nad83csrs.SHIFT_FILE = script_base + '/share/NAD83v6VG.tif'
nad83csrs.ITRF_FILE = script_base + '/share/itrf.csv'

pt = loads(geom)
t = nad83csrs.Transformer(ffrom, efrom, eto, fromsrs, tosrs)
tries = 5;
exc = None
while tries >= 0:
	try:
		x, y, z, bounds = t.transformPoints([pt.x], [pt.y], [pt.z])
		return dumps(Point(x[0], y[0], z[0]))
	except Exception, e:
		exc = e
	tries -= 1
if exc:
	raise exc

$$ LANGUAGE 'plpythonu';

CREATE OR REPLACE FUNCTION ToNAD83CSRS(geom geometry, ffrom text, efrom float8, eto float8, fromsrid integer, tosrid integer)
RETURNS geometry
AS $$
BEGIN
	return st_setsrid(st_geomfromewkb(_ToNAD83CSRS(st_asewkb(geom), ffrom, efrom, eto, 
		(select proj4text from spatial_ref_sys where srid=fromsrid), 
		(select proj4text from spatial_ref_sys where srid=tosrid)
	)), tosrid);
END;
$$ LANGUAGE 'plpgsql';

-- Transforms a PostGIS geometry in any coordinate reference system at any epoch and reference frame
-- to NAD83(CSRS) at any epoch.
-- 
-- geoms 	An array of 3D point geometries in WKB form.
-- ffrom	The origin reference frame. This will be something like 'ITRF90'.
-- efrom	The origin epoch. Epochs are expressed as decimal years.
-- eto		The destination epoch.
-- fromSRS  The proj string of the origin geometry.
-- toSRS    The proj string of the destination geometry.
CREATE OR REPLACE FUNCTION ToNAD83CSRS(geoms bytea[], ffrom text, efrom float8, eto float8, fromsrs text, tosrs text)
RETURNS bytea[]
AS $$

global fromsrs

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

count = len(geoms)
x = [None] * count
y = [None] * count
z = [None] * count

for i in range(count):
	pt = loads(geoms[i])
	x[i] = pt.x
	y[i] = pt.y
	z[i] = pt.z

t = nad83csrs.Transformer(ffrom, efrom, eto, fromsrs, tosrs)
x, y, z, bounds = t.transformPoints(x, y, z)

result = [None] * count
for i in range(count):
	result[i] = dumps(Point(x[i], y[i], z[i]))

return result

$$ LANGUAGE 'plpythonu';
