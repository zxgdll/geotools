cmake_minimum_required(VERSION 3.0)

if (${APPLE})
	set (CMAKE_C_COMPILER /usr/local/bin/gcc-6 CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/g++-6 CACHE STRING "C++ compiler" FORCE)
endif()

enable_language(C)
enable_language(CXX)

find_package(OpenMP)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_CXX_STANDARD 11)

project (geotools)

link_directories(/usr/local/lib /usr/lib)
include_directories(include /usr/local/include /usr/include /usr/local/include/gdal /usr/local/include/geos /usr/include/gdal /usr/include/geos)


add_executable (flood src/flood.cpp)
target_link_libraries (flood gdal)



add_executable (rasterstats src/rasterstats.cpp)
target_link_libraries (rasterstats gdal)

add_executable (mosaic src/mosaic.cpp)
target_link_libraries (mosaic gdal boost_system boost_filesystem)



add_library (treetops SHARED src/treetops.cpp)
target_link_libraries (treetops gdal)
add_executable (treetops-app src/apps/treetops.cpp)
target_link_libraries (treetops-app treetops gdal)



add_library(lasgrid SHARED src/lasgrid.cpp)
target_link_libraries (lasgrid las gdal las laszip boost_system boost_filesystem)
add_executable (lasgrid-app src/apps/lasgrid.cpp)
target_link_libraries (lasgrid-app PUBLIC lasgrid)

add_executable (lasboundary src/lasboundary.cpp)
target_link_libraries (lasboundary PUBLIC las gdal las laszip geos boost_system boost_filesystem)

add_library(lasvalidate SHARED src/lasvalidate.cpp)
target_link_libraries (lasvalidate las laszip geos)
add_executable (lasvalidate-app src/apps/lasvalidate.cpp)
target_link_libraries (lasvalidate-app lasvalidate)

add_executable (lasstats src/lasstats.cpp)
target_link_libraries (lasstats gdal las laszip geos boost_system boost_filesystem)

add_executable (lasmerge src/lasmerge.cpp)
target_link_libraries (lasmerge las laszip boost_system boost_filesystem)

add_executable (lasreclass src/lasreclass.cpp)
target_link_libraries (lasreclass las laszip boost_system boost_filesystem)

add_executable (lasclip src/lasclip.cpp)
target_link_libraries (lasclip gdal las laszip geos boost_system boost_filesystem)

add_executable (las2csrs src/las2csrs.cpp)
target_link_libraries (las2csrs gdal las laszip proj boost_system boost_filesystem)



install (TARGETS flood DESTINATION bin)

install (TARGETS treetops DESTINATION lib)
install (TARGETS treetops-app DESTINATION bin)
install (FILES include/treetops.hpp DESTINATION include)

install (TARGETS rasterstats DESTINATION bin)

install (TARGETS mosaic DESTINATION bin)

install (TARGETS lasgrid DESTINATION lib)
install (TARGETS lasgrid-app DESTINATION bin)
install (FILES include/lasgrid.hpp DESTINATION include)

install (TARGETS lasvalidate DESTINATION lib)
install (TARGETS lasvalidate-app DESTINATION bin)
install (FILES include/lasvalidate.hpp DESTINATION include)

install (TARGETS lasstats DESTINATION bin)

install (TARGETS lasmerge DESTINATION bin)

install (TARGETS lasreclass DESTINATION bin)

install (TARGETS lasclip DESTINATION bin)

install (TARGETS las2csrs DESTINATION bin)
install (DIRECTORY share/ DESTINATION share PATTERN "share/*")

