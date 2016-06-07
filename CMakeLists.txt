cmake_minimum_required (VERSION 3.1.0)
set (CMAKE_CXX_STANDARD 11)

project (geotools)

link_directories(/usr/local/lib)

add_executable (lasgrid src/lasgrid.cpp)

target_include_directories (lasgrid PUBLIC include)
target_include_directories (lasgrid PUBLIC /usr/local/include)
target_link_libraries (lasgrid gdal las boost_system boost_filesystem)

add_executable (lasvalidate src/lasvalidate.cpp)

target_include_directories (lasvalidate PUBLIC include)
target_include_directories (lasvalidate PUBLIC /usr/local/include)
target_link_libraries (lasvalidate las geos)

add_executable (las2csrs src/las2csrs.cpp)

target_include_directories (las2csrs PUBLIC include)
target_include_directories (las2csrs PUBLIC /usr/local/include)
target_link_libraries (las2csrs gdal las proj boost_system boost_filesystem)

add_executable (lasstats src/lasstats.cpp)

target_include_directories (lasstats PUBLIC include)
target_include_directories (lasstats PUBLIC /usr/local/include)
target_link_libraries (lasstats gdal las geos boost_system boost_filesystem)

add_executable (lasmerge src/lasmerge.cpp)

target_include_directories (lasmerge PUBLIC include)
target_include_directories (lasmerge PUBLIC /usr/local/include)
target_link_libraries (lasmerge las boost_system boost_filesystem)

add_executable (lasreclass src/lasreclass.cpp)

target_include_directories (lasreclass PUBLIC include)
target_include_directories (lasreclass PUBLIC /usr/local/include)
target_link_libraries (lasreclass las boost_system boost_filesystem)

add_executable (lasclip src/lasclip.cpp)

target_include_directories (lasclip PUBLIC include)
target_include_directories (lasclip PUBLIC /usr/local/include)
target_link_libraries (lasclip gdal las geos boost_system boost_filesystem)

IF (OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
ENDIF ()

add_executable (mosaic src/mosaic.cpp)

target_include_directories (mosaic PUBLIC include)
target_include_directories (mosaic PUBLIC /usr/local/include)
target_link_libraries (mosaic gdal boost_system boost_filesystem)

install (TARGETS lasgrid DESTINATION bin)
install (TARGETS lasvalidate DESTINATION bin)
install (TARGETS las2csrs DESTINATION bin)
install (DIRECTORY share/ DESTINATION share PATTERN "share/*")
install (TARGETS lasstats DESTINATION bin)
install (TARGETS lasmerge DESTINATION bin)
install (TARGETS lasreclass DESTINATION bin)
install (TARGETS lasclip DESTINATION bin)
install (TARGETS mosaic DESTINATION bin)
