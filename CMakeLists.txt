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
include_directories(include /usr/local/include /usr/include /usr/include/gdal)

add_executable (flood src/flood.cpp)
target_link_libraries (flood gdal)

add_executable (rasterstats src/rasterstats.cpp)
target_link_libraries (rasterstats gdal)

add_executable (treetops src/treetops.cpp)
target_link_libraries (treetops gdal)

add_executable (lasgrid src/lasgrid.cpp)
target_link_libraries (lasgrid PUBLIC las gdal las laszip boost_system boost_filesystem)

#add_library(lasvalidate SHARED src/lasvalidate.cpp)
#target_link_libraries (lasvalidate las laszip geos)

#add_executable (lasvalidate-app src/apps/lasvalidate.cpp)
#target_link_libraries (lasvalidate-app lasvalidate las geos)

add_executable (las2csrs src/las2csrs.cpp)
target_link_libraries (las2csrs gdal las laszip proj boost_system boost_filesystem)

add_executable (lasstats src/lasstats.cpp)
target_link_libraries (lasstats gdal las laszip geos boost_system boost_filesystem)

add_executable (lasmerge src/lasmerge.cpp)
target_link_libraries (lasmerge las laszip boost_system boost_filesystem)

add_executable (lasreclass src/lasreclass.cpp)
target_link_libraries (lasreclass las laszip boost_system boost_filesystem)

add_executable (lasclip src/lasclip.cpp)
target_link_libraries (lasclip gdal las laszip geos boost_system boost_filesystem)

add_executable (mosaic src/mosaic.cpp)
target_link_libraries (mosaic gdal boost_system boost_filesystem)

install (TARGETS lasgrid DESTINATION bin)
#install (TARGETS lasvalidate DESTINATION lib)
#install (TARGETS lasvalidate-app DESTINATION bin)
install (FILES include/lasvalidate.hpp DESTINATION include)
install (TARGETS las2csrs DESTINATION bin)
install (DIRECTORY share/ DESTINATION share PATTERN "share/*")
install (TARGETS lasstats DESTINATION bin)
install (TARGETS lasmerge DESTINATION bin)
install (TARGETS lasreclass DESTINATION bin)
install (TARGETS lasclip DESTINATION bin)
install (TARGETS mosaic DESTINATION bin)
install (TARGETS rasterstats DESTINATION bin)
install (TARGETS flood DESTINATION bin)
install (TARGETS treetops DESTINATION bin)

