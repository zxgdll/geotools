DELETE FROM spatial_ref_sys WHERE srid>99990;

-- SRID for canadian peatland atlas
-- Geological Survey of Canada Open File 3834
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 200000, 'epsg', 200000, '+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +ellps=clrk66 +units=m +no_defs','');

-- 26912 + 5713
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 99991, 'epsg', 26912, '+init=EPSG:26912 +geoidgrids=HT2_0.gtx', 
	'PROJCS["NAD83 / UTM zone 12N",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","26912"],AXIS["Easting",EAST],AXIS["Northing",NORTH]],
	VERT_CS["Canadian Geodetic Vertical Datum of 1928",VERT_DATUM["Canadian Geodetic Vertical Datum of 1928",2005,AUTHORITY["EPSG","5114"]],UNIT["m",1.0],AXIS["Gravity-related height",UP],AUTHORITY["EPSG","5713"]]');

-- 26712 + 5713
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100000, 'epsg', 26712, '+init=EPSG:26712 +geoidgrids=HT2_0.gtx', 
	'PROJCS["NAD27 / UTM zone 12N",GEOGCS["NAD27",DATUM["North_American_Datum_1927",SPHEROID["Clarke 1866",6378206.4,294.9786982139006,AUTHORITY["EPSG","7008"]],AUTHORITY["EPSG","6267"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4267"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","26712"]],
	VERT_CS["Canadian Geodetic Vertical Datum of 1928",VERT_DATUM["Canadian Geodetic Vertical Datum of 1928",2005,AUTHORITY["EPSG","5114"]],UNIT["m",1.0],AXIS["Gravity-related height",UP],AUTHORITY["EPSG","5713"]]');

-- 26712 + CGG2010
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100001, 'epsg', 26712, '+init=EPSG:26712 +geoidgrids=CGG2010n83.gtx', 
	'PROJCS["NAD27 / UTM zone 12N",GEOGCS["NAD27",DATUM["North_American_Datum_1927",SPHEROID["Clarke 1866",6378206.4,294.9786982139006,AUTHORITY["EPSG","7008"]],AUTHORITY["EPSG","6267"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4267"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","26712"]]');

-- 26912 + CGG2000
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100002, 'epsg', 26912, '+init=EPSG:26912 +geoidgrids=CGG2000n83.gtx', 
	'PROJCS["NAD83 / UTM zone 12N",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","26912"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 26912 + NGSD95n83
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 99999, 'epsg', 26912, '+init=EPSG:26912 +geoidgrids=NGSD95n83.gtx', 
	'PROJCS["NAD83 / UTM zone 12N",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","26912"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 26912 + GSD95
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100020, 'epsg', 26912, '+init=EPSG:26912 +geoidgrids=GSD95.gtx', 
	'PROJCS["NAD83 / UTM zone 12N",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","26912"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 26912 + NGSD91n83
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100006, 'epsg', 26912, '+init=EPSG:26912 +geoidgrids=NGSD91n83.gtx', 
	'PROJCS["NAD83 / UTM zone 12N",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","26912"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 26912 + GSD91
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100022, 'epsg', 26912, '+init=EPSG:26912 +geoidgrids=GSD91.gtx', 
	'PROJCS["NAD83 / UTM zone 12N",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","26912"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 4326 + NGSD91n83
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100011, 'epsg', 4326, '+init=4326 +geoidgrids=NGSD91n83.gtx', 
	'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]');

-- 32612 + NGSD95n83
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100005, 'epsg', 32612, '+init=EPSG:32612 +geoidgrids=NGSD95n83.gtx', 
	'PROJCS["WGS 84 / UTM zone 12N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32612"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 32612 + GSD95
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100021, 'epsg', 32612, '+init=EPSG:32612 +geoidgrids=GSD95.gtx', 
	'PROJCS["WGS 84 / UTM zone 12N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32612"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 32612 + 5713
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 99992, 'epsg', 32612, '+init=EPSG:32612 +geoidgrids=HT2_0.gtx', 
	'PROJCS["WGS 84 / UTM zone 12N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32612"],AXIS["Easting",EAST],AXIS["Northing",NORTH]],VERT_CS["Canadian Geodetic Vertical Datum of 1928",VERT_DATUM["Canadian Geodetic Vertical Datum of 1928",2005,AUTHORITY["EPSG","5114"]],UNIT["m",1.0],AXIS["Gravity-related height",UP],AUTHORITY["EPSG","5713"]]');

-- 32612 + 6647
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 99993, 'epsg', 32612, '+init=EPSG:32612 +geoidgrids=CGG2013n83.gtx', 
	'PROJCS["WGS 84 / UTM zone 12N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32612"],AXIS["Easting",EAST],AXIS["Northing",NORTH]],
	VERT_CS["Canadian Geodetic Vertical Datum of 2013",VERT_DATUM["Canadian Geodetic Vertical Datum of 2013",2005,AUTHORITY["EPSG","1127"]],UNIT["m",1.0],AXIS["Gravity-related height",UP],AUTHORITY["EPSG","6647"]]');

-- 2956 + 6647 (also 6655)
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100007, 'epsg', 2956, '+init=EPSG:2956 +geoidgrids=CGG2013n83.gtx', 
	'PROJCS["NAD83(CSRS) / UTM zone 12N",GEOGCS["NAD83(CSRS)",DATUM["NAD83_Canadian_Spatial_Reference_System",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6140"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4617"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","2956"],AXIS["Easting",EAST],AXIS["Northing",NORTH]],
	VERT_CS["Canadian Geodetic Vertical Datum of 2013",VERT_DATUM["Canadian Geodetic Vertical Datum of 2013",2005,AUTHORITY["EPSG","1127"]],UNIT["m",1.0],AXIS["Gravity-related height",UP],AUTHORITY["EPSG","6647"]]');

-- 2956 + 5713
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100008, 'epsg', 2956, '+init=EPSG:2956 +geoidgrids=HT2_0.gtx', 
	'PROJCS["NAD83(CSRS) / UTM zone 12N",GEOGCS["NAD83(CSRS)",DATUM["NAD83_Canadian_Spatial_Reference_System",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6140"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4617"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","2956"],AXIS["Easting",EAST],AXIS["Northing",NORTH]],
	VERT_CS["Canadian Geodetic Vertical Datum of 1928",VERT_DATUM["Canadian Geodetic Vertical Datum of 1928",2005,AUTHORITY["EPSG","5114"]],UNIT["m",1.0],AXIS["Gravity-related height",UP],AUTHORITY["EPSG","5713"]]');

-- 2956 + NGSD95n83
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100012, 'epsg', 2956, '+init=EPSG:2956 +geoidgrids=NGSD95n83.gtx', 
	'PROJCS["NAD83(CSRS) / UTM zone 12N",GEOGCS["NAD83(CSRS)",DATUM["NAD83_Canadian_Spatial_Reference_System",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6140"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4617"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","2956"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 2956 + NGSD91n83
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100014, 'epsg', 2956, '+init=EPSG:2956 +geoidgrids=NGSD91n83.gtx', 
	'PROJCS["NAD83(CSRS) / UTM zone 12N",GEOGCS["NAD83(CSRS)",DATUM["NAD83_Canadian_Spatial_Reference_System",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6140"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4617"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","2956"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 2956 + cgg2000
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100013, 'epsg', 2956, '+init=EPSG:2956 +geoidgrids=CGG2000n83.gtx', 
	'PROJCS["NAD83(CSRS) / UTM zone 12N",GEOGCS["NAD83(CSRS)",DATUM["NAD83_Canadian_Spatial_Reference_System",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6140"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4617"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","2956"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 4326 + 5713
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 99994, 'epsg', 4326, '+init=EPSG:4326 +geoidgrids=HT2_0.gtx', 
	'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]], 
	VERT_CS["Canadian Geodetic Vertical Datum of 1928",VERT_DATUM["Canadian Geodetic Vertical Datum of 1928",2005,AUTHORITY["EPSG","5114"]],UNIT["m",1.0],AXIS["Gravity-related height",UP],AUTHORITY["EPSG","5713"]]');

-- 4326 + 6647
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100004, 'epsg', 4326, '+init=EPSG:4326 +geoidgrids=CGG2013n83.gtx', 
	'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]], 
	VERT_CS["Canadian Geodetic Vertical Datum of 2013",VERT_DATUM["Canadian Geodetic Vertical Datum of 2013",2005,AUTHORITY["EPSG","1127"]],UNIT["m",1.0],AXIS["Gravity-related height",UP],AUTHORITY["EPSG","6647"]]');

-- 2956 + CGG2010 
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 99996, 'epsg', 2956, '+init=EPSG:2956 +geoidgrids=CGG2010n83.gtx', 
	'PROJCS["NAD83(CSRS) / UTM zone 12N",GEOGCS["NAD83(CSRS)",DATUM["NAD83_Canadian_Spatial_Reference_System",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6140"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4617"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","2956"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 32612 + CGG2010 
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100009, 'epsg', 32612, '+init=EPSG:32612 +geoidgrids=CGG2010n83.gtx', 
	'PROJCS["NAD83(CSRS) / UTM zone 12N",GEOGCS["NAD83(CSRS)",DATUM["NAD83_Canadian_Spatial_Reference_System",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6140"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4617"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","2956"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 32612 + CGG2000 
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 100010, 'epsg', 32612, '+init=EPSG:32612 +geoidgrids=CGG2000n83.gtx', 
	'PROJCS["NAD83(CSRS) / UTM zone 12N",GEOGCS["NAD83(CSRS)",DATUM["NAD83_Canadian_Spatial_Reference_System",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6140"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4617"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","2956"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]');

-- 2956 + CGVD28
INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 99997, 'epsg', 2956, '+init=EPSG:2956 +geoidgrids=HT2_0.gtx', 
	'PROJCS["NAD83(CSRS) / UTM zone 12N",GEOGCS["NAD83(CSRS)",DATUM["NAD83_Canadian_Spatial_Reference_System",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6140"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4617"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","2956"],AXIS["Easting",EAST],AXIS["Northing",NORTH]],
	VERT_CS["Canadian Geodetic Vertical Datum of 1928",VERT_DATUM["Canadian Geodetic Vertical Datum of 1928",2005,AUTHORITY["EPSG","5114"]],UNIT["m",1.0],AXIS["Gravity-related height",UP],AUTHORITY["EPSG","5713"]]');


-- -0.991,1.9072,0.5129,-0.0257899075194932,-0.0096500989602704,-0.0116599432323421,0,1
-- ,VERT_CS["Canadian Geodetic Vertical Datum of 2010",VERT_DATUM["Canadian Geodetic Vertical Datum of 2010",2005,AUTHORITY["EPSG","1127"]],UNIT["m",1.0],AXIS["Gravity-related height",UP],AUTHORITY["EPSG","99995"]]

--old  | POINT Z (493377.131 6514034.137 225.059002491229) | POINT Z (493377.131 6514034.137 225.01) | 
--new  | POINT Z (493377.2 6514034.145 253.444903905764)   |                                         | POINT Z (493377.2 6514034.145 224.285)

--POINT Z (493377.131 6514034.137 225.01) --> POINT Z (493377.131 6514034.137 225.363034643733)
--POINT Z (493377.2 6514034.145 224.285) --> POINT Z (493377.2 6514034.145 224.158305729835)
