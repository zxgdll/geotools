 CREATE EXTENSION IF NOT EXISTS plpythonu;


-- Transforms a PostGIS geometry in any coordinate reference system at any epoch and reference frame
-- to NAD83(CSRS) at any epoch.
-- 
-- geom 	A 3D point geometry in WKB form.
-- ffrom	The origin reference frame. This will be something like 'ITRF90'.
-- efrom	The origin epoch. Epochs are expressed as decimal years.
-- eto		The destination epoch.
CREATE OR REPLACE FUNCTION ToNAD83CSRS(geom bytea, ffrom text, efrom double precision, eto double precision, type text, zone integer)
RETURNS bytea
AS $$

from nad83csrs import transform
from shapely.geometry import Point
from shapely.wkb import loads, dumps

pt = loads(geom)
x, y, z = transform(pt.x, pt.y, pt.z, ffrom, efrom, eto, type, zone)
return dumps(Point(x, y, z))

$$ LANGUAGE 'plpythonu';


-- This a temporary method for configuring the init.
CREATE OR REPLACE FUNCTION ToNAD83CSRS_Init()
RETURNS VOID
AS $$

SELECT ToNAD83CSRS_Init('/Users/robskelly/Documents/geotools/scripts', 
	'/Users/robskelly/Documents/geotools/scripts/NAD83v6VG.tif', 
	'/Users/robskelly/Documents/geotools/scripts/itrf.csv')

$$ LANGUAGE 'sql';


-- Initializes the NAD83(CSRS) transformation package, which requires loading some 
-- external data and initializing some paths. If the _Cleanup function is not
-- called subsequently, a memory leak will occur!
CREATE OR REPLACE FUNCTION ToNAD83CSRS_Init(script_path text, shift_file text, itrf_file text)
RETURNS VOID
AS $$

import sys

# Set the path for the script.
sys.path.append(script_path)

import nad83csrs

# Give the module a reference to PostgreSQL's persistent storage.
nad83csrs.DATA = SD

# Set the paths for required data files.
nad83csrs.SHIFT_FILE = shift_file
nad83csrs.ITRF_FILE = itrf_file

nad83csrs.init()

$$ LANGUAGE 'plpythonu';


-- Cleans up resources used by the NAD83(CSRS) transformation routines.
-- This must be called afterwards or a memory leak will occur!
CREATE OR REPLACE FUNCTION ToNAD83CSRS_Cleanup()
RETURNS VOID
AS $$

import nad83csrs

nad83csrs.cleanup()

$$ LANGUAGE 'plpythonu';
