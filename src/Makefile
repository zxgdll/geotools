LIB=-L/opt/lib
CXXFLAGS=-I../include -I/opt/include
CFLAGS=-I../include -I/opt/include
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
	
proj_ellips:
	$(CC) $(CFLAGS) -Wall -g proj_ellips.c $(LIB) -o ../bin/proj_ellips

all: mosaic, lasstats, lasmerge, lasgrid
debug: mosaic_debug, laststats_debug, lasmerge_debug, lasgrid_debug

clean:
	rm ../bin/*
