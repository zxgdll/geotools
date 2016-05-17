CXXFLAGS=-std=c++11 -Wall -I../include -I/usr/local/include -I/usr/include 

# To deal with variable path prefixes
CXXFLAGS += -I/usr/local/include/gdal -I/usr/include/gdal 
CXXFLAGS += -I/usr/include/local/liblas -I/usr/include/liblas 
CXXFLAGS += -I/usr/include/local/eigen3 -I/usr/include/eigen3

# Using a library (like GEOS trunk) that is installed in /opt
#CXXFLAGS := -I/opt/include $(CXXFLAGS)

LIB=-L/usr/local/lib -L/usr/lib
CC=gcc
CPP=g++

.PHONY: interp

all: mosaic lasmerge lasgrid las2csrs lasclip interp lasstats

debug: CXXFLAGS += -DDEBUG -g
debug: CFLAGS += -DDEBUG -g
debug: $(target)

release: CXXFLAGS += -O3
release: CFLAGS += -O3
release: $(target)

mosaic:
	$(CPP) $(CXXFLAGS) -fopenmp mosaic.cpp $(LIB) -lgdal -o ../bin/mosaic

treecrowns:
	$(CPP) $(CXXFLAGS) treecrowns.cpp $(LIB) -lgdal -o ../bin/treecrowns

lasstats:
	$(CPP) $(CXXFLAGS) lasstats.cpp $(LIB) -llas -lgdal -lproj -lgeos -o ../bin/lasstats

lasmerge:
	$(CPP) $(CXXFLAGS) lasmerge.cpp $(LIB) -lboost_system -llas -o ../bin/lasmerge

lasgrid:
	$(CPP) $(CXXFLAGS) lasgrid.cpp $(LIB) -llas -lgdal -lboost_system -o ../bin/lasgrid

las2csrs:
	$(CPP) $(CXXFLAGS) las2csrs.cpp $(LIB) -llas -lgdal -lboost_system -lboost_filesystem -lproj -o ../bin/las2csrs

lasclip:
	$(CPP) $(CXXFLAGS) lasclip.cpp $(LIB) -llas -lgeos -lgdal -lboost_system -o ../bin/lasclip

#rastfit:
#	$(CPP) $(CXXFLAGS) -DWITH_CGAL interp/Interpolators.cpp rastfit.cpp $(LIB) -lgdal -lmpfr -lgmp -lCGAL -lCGAL_Core -lboost_thread-mt -o ../bin/rastfit

interp:
	$(CPP) $(CXXFLAGS) -fopenmp interp/Interpolators.cpp interp.cpp $(LIB) -lgomp -lgdal -lmpfr -lgmp -lCGAL -lCGAL_Core -o ../bin/interp

lasvalidate:
	$(CPP) $(CXXFLAGS) lasvalidate.cpp $(LIB) -llas -lgeos -o ../bin/lasvalidate

lasboundary:
	$(CPP) $(CXXFLAGS) -fopenmp lasboundary.cpp $(LIB) -llas -lgomp -lmpfr -lgmp -lgdal -lgeos -lboost_system -lboost_filesystem -o ../bin/lasboundary

#proj_ellips:
#	$(CC) $(CFLAGS) proj_ellips.c $(LIB) -o ../bin/proj_ellips

#findmin:
#	$(CPP) $(CXXFLAGS) utils/findmin.cpp $(LIB) -lgdal -o ../bin/findmin

#srtm_lakes:
#	$(CPP) $(CXXFLAGS) utils/srtm_lakes.cpp $(LIB) -lgdal -o ../bin/srtm_lakes

clean:
	rm ../bin/*
