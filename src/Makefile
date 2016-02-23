LIB=-L/usr/lib -L/usr/lib/qt4
CXXFLAGS=-I../include -I/usr/include -I/usr/include/gdal -I/usr/include/qt4 -I/usr/include/qt4/Qt -I/usr/include/qwt -std=c++11
CFLAGS=-I../include -I/usr/include -std=c++11
CPP=g++
CC=gcc

.PHONY: interp

all: mosaic lasstats lasmerge lasgrid las2csrs rastfit lasclip lasboundary lasreclass
debug: mosaic_debug laststats_debug lasmerge_debug lasgrid_debug las2csrs_debug rastfit_debug lasclip_debug lasboundary_debug lasreclass_debug

clean:
	rm ../bin/*

mosaic:
	$(CPP) $(CXXFLAGS) -Wall -O3 mosaic.cpp $(LIB) -fopenmp -lgdal -o ../bin/mosaic

mosaic_debug:
	$(CPP) $(CXXFLAGS) -Wall -g mosaic.cpp $(LIB) -fopenmp -lgdal -o ../bin/mosaic
	
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
	
lasboundary:
	$(CPP) $(CXXFLAGS) -Wall -O3 lasboundary.cpp $(LIB) -llas -lgdal -lboost_system -lboost_filesystem -o ../bin/lasboundary

lasboundary_debug:
	$(CPP) $(CXXFLAGS) -Wall -g lasboundary.cpp $(LIB) -llas -lgdal -lboost_system -lboost_filesystem -o ../bin/lasboundary
	
lasclip:
	$(CPP) $(CXXFLAGS) -Wall -O3 lasclip.cpp $(LIB) -llas -lgeos -lgdal -lboost_system -o ../bin/lasclip

lasclip_debug:
	$(CPP) $(CXXFLAGS) -Wall -g lasclip.cpp $(LIB) -llas -lgeos -lgdal -lboost_system -o ../bin/lasclip
	
lasreclass:
	$(CPP) $(CXXFLAGS) -Wall -O3 lasreclass.cpp $(LIB) -llas -lboost_system -o ../bin/lasreclass

lasreclass_debug:
	$(CPP) $(CXXFLAGS) -Wall -g lasreclass.cpp $(LIB) -llas -lboost_system -o ../bin/lasreclass
	
las2csrs:
	$(CPP) $(CXXFLAGS) -Wall -O3 las2csrs.cpp $(LIB) -llas -lgdal -lboost_system -lboost_filesystem -lproj -o ../bin/las2csrs

las2csrs_debug:
	$(CPP) $(CXXFLAGS) -Wall -g las2csrs.cpp $(LIB) -llas -lgdal -lboost_system -lboost_filesystem -lproj -o ../bin/las2csrs
	
rastfit:
	$(CPP) $(CXXFLAGS) -Wall -O3 interp/Interpolators.cpp rastfit.cpp $(LIB) -lQtCore -lQtGui -lqwt -lgdal -lmpfr -lgmp -lCGAL -lCGAL_Core -lboost_thread -o ../bin/rastfit

rastfit_debug:
	$(CPP) $(CXXFLAGS) -Wall -g interp/Interpolators.cpp rastfit.cpp $(LIB) -lgdal -lmpfr -lgmp -lCGAL -lCGAL_Core -lboost_thread -o ../bin/rastfit
	
interp:
	$(CPP) $(CXXFLAGS) -Wall -O3 interp/Interpolators.cpp interp.cpp $(LIB) -lQtCore -lQtGui -lqwt -lgdal -lmpfr -lgmp -lCGAL -lCGAL_Core -lboost_thread -o ../bin/interp

interp_debug:
	$(CPP) $(CXXFLAGS) -Wall -g interp/Interpolators.cpp interp.cpp $(LIB) -lQtCore -lQtGui -lqwt -lgdal -lmpfr -lgmp -lCGAL -lCGAL_Core -lboost_thread -o ../bin/interp
	
proj_ellips:
	$(CC) $(CFLAGS) -Wall -g proj_ellips.c $(LIB) -o ../bin/proj_ellips

