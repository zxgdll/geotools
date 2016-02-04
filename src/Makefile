LIB=-L/usr/lib
CXXFLAGS=-I../include -I/usr/include -I/usr/include/gdal -std=c++11
CFLAGS=-I../include -I/usr/include -std=c++11
CPP=g++
CC=gcc

# -lboost_system -lboost_iostreams -lboost_program_options -lboost_thread -lproj -lboost_filesystem -lgdal -lgeos -llas -lgmp -lmpfr

mosaic:
	$(CPP) $(CXXFLAGS) -Wall -O3 mosaic.cpp $(LIB) -lgdal -o ../bin/mosaic

mosaic_debug:
	$(CPP) $(CXXFLAGS) -Wall -g mosaic.cpp $(LIB) -lgdal -o ../bin/mosaic
	
lasstats:
	$(CPP) $(CXXFLAGS) -Wall -O3 lasstats.cpp $(LIB) -llas -lgdal -lproj -lgeos -o ../bin/lasstats

lasstats_debug:
	$(CPP) $(CXXFLAGS) -Wall -g lasstats.cpp $(LIB) -llas -lgdal -lproj -lgeos -o ../bin/lasstats
	
lasmerge:
	$(CPP) $(CXXFLAGS) -Wall -O3 lasmerge.cpp $(LIB) -llas -lgdal -lgeos -lboost_system -o ../bin/lasmerge

lasmerge_debug:
	$(CPP) $(CXXFLAGS) -Wall -g lasmerge.cpp $(LIB) -llas -lgdal -lgeos  -o ../bin/lasmerge
	
lasgrid:
	$(CPP) $(CXXFLAGS) -Wall -O3 lasgrid.cpp $(LIB) -llas -lgdal -lboost_system -o ../bin/lasgrid

lasgrid_debug:
	$(CPP) $(CXXFLAGS) -Wall -g lasgrid.cpp $(LIB) -llas -lgdal -lboost_system -o ../bin/lasgrid
	
las2csrs:
	$(CPP) $(CXXFLAGS) -Wall -O3 las2csrs.cpp $(LIB) -llas -lgdal -lboost_system -lboost_filesystem -lproj -o ../bin/las2csrs

las2csrs_debug:
	$(CPP) $(CXXFLAGS) -Wall -g las2csrs.cpp $(LIB) -llas -lgdal -lboost_system -lboost_filesystem -lproj -o ../bin/las2csrs
	
rastfit:
	$(CPP) $(CXXFLAGS) -Wall -O3 rastfit.cpp $(LIB) -lgdal -lgmp -lCGAL -lCGAL_Core -lboost_thread -o ../bin/rastfit

rastfit_debug:
	$(CPP) $(CXXFLAGS) -Wall -g rastfit.cpp $(LIB) -lgdal -lgmp -lCGAL -lCGAL_Core -lboost_thread -o ../bin/rastfit
	
proj_ellips:
	$(CC) $(CFLAGS) -Wall -g proj_ellips.c $(LIB) -o ../bin/proj_ellips

all: mosaic, lasstats, lasmerge, lasgrid
debug: mosaic_debug, laststats_debug, lasmerge_debug, lasgrid_debug

clean:
	rm ../bin/*
