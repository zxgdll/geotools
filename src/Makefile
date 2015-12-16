LIB=-L/opt/geos/lib -L/opt/gdal/lib -L/opt/liblas/lib -lboost_system -lboost_iostreams -lboost_program_options -lboost_thread -lproj -lboost_filesystem -lgdal -llas -lgmp -lmpfr
CXXFLAGS=-I../include -I/opt/gdal/include -I/opt/liblas/include -I/opt/geos/include
CFLAGS=-I../include -I/opt/gdal/include -I/opt/liblas/include -I/opt/geos/include
CPP=g++
CC=gcc

all: mosaic
debug: mosaic_debug

mosaic:
	$(CPP) $(CXXFLAGS) -Wall -O3 mosaic.cpp $(LIB) -o ../bin/mosaic

mosaic_debug:
	$(CPP) $(CXXFLAGS) -Wall -g mosaic.cpp $(LIB) -o ../bin/mosaic
	
proj_ellips:
	$(CC) $(CFLAGS) -Wall -g proj_ellips.c $(LIB) -o ../bin/proj_ellips

clean:
	rm ../bin/*
