cmake_minimum_required(VERSION 3.0)

set (WITH_GUI false CACHE BOOL "If true, will build the Qt GUI.")

if(${WITH_GUI})
	add_definitions(-DWITH_GUI)
endif()

if (${APPLE})
	if(${WITH_GUI})
		set (Qt5_DIR /usr/local/Cellar/qt5/5.6.1-1/lib/cmake/Qt5)
	endif()
	set (CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE STRING "C++ compiler" FORCE)
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
endif()

enable_language(C)
enable_language(CXX)

find_package(OpenMP)
if(${WITH_GUI})
	find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets)
endif()

set (CMAKE_AUTOUIC ON)
set (CMAKE_AUTOMOC ON)
set (CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_AUTOMOC_OPTIONS "-Iinclude")

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_CXX_STANDARD 11)

set (CPACK_PACKAGE_NAME "GeoTools")
set (CPACK_PACKAGE_VENDOR "dijital.ca")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "GeoTools - Tools for working with geodata.")
set (CPACK_PACKAGE_VERSION "1.0.0")
set (CPACK_PACKAGE_VERSION_MAJOR "1")
set (CPACK_PACKAGE_VERSION_MINOR "0")
set (CPACK_PACKAGE_VERSION_PATCH "0")
set (CPACK_PACKAGE_INSTALL_DIRECTORY "geotools")

set (CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "GeoTools Applications")
set (CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "GeoTools Libraries")
set (CPACK_COMPONENT_HEADERS_DISPLAY_NAME "GeoTools Header Files")

set (CPACK_COMPONENT_APPLICATIONS_DESCRIPTION "GeoTools applications.")
set (CPACK_COMPONENT_LIBRARIES_DESCRIPTION "Libraries for extending GeoTools.")
set (CPACK_COMPONENT_HEADERS_DESCRIPTION "C/C++ header files for use with GeoTools.")

set (CPACK_COMPONENT_LIBRARIES_INSTALL_TYPES Full)
set (CPACK_COMPONENT_APPLICATIONS_INSTALL_TYPES Full User)
set (CPACK_COMPONENT_HEADERS_INSTALL_TYPES Full)

set (CPACK_ALL_INSTALL_TYPES Full User)
set (CPACK_COMPONENTS_ALL applications libraries headers)
set (CPACK_INSTALL_CMAKE_PROJECTS /Users/robskelly/Documents/geotools/makefiles/pointstats-app;pointstats-app;ALL;/)

include (CPack)

project (geotools)

link_directories(/usr/local/opt/sqlite/lib /usr/local/lib /usr/lib)
include_directories(/usr/local/opt/sqlite/include include /usr/local/include /usr/include /usr/include/x86_64-linux-gnu/qt5 /usr/local/include/gdal /usr/local/include/geos /usr/include/gdal /usr/include/geos)


add_library (geotools SHARED src/geotools.cpp)
add_library (util SHARED src/util.cpp)
add_library (raster SHARED src/raster.cpp)
target_link_libraries (raster gdal util geotools)


add_executable (flood src/flood.cpp)
target_link_libraries (flood gdal raster util geotools)


add_executable (rasterstats src/rasterstats.cpp)
target_link_libraries (rasterstats gdal raster util geotools)


add_library (mosaic SHARED src/mosaic.cpp)
target_link_libraries (mosaic gdal boost_system boost_filesystem raster util geotools)
if(${WITH_GUI})
	add_executable (mosaic-app src/apps/mosaic.cpp src/apps/mosaic_ui.cpp)
	target_link_libraries (mosaic-app PUBLIC mosaic Qt5::Gui Qt5::Core Qt5::Widgets)
else()
	add_executable (mosaic-app src/apps/mosaic.cpp)
	target_link_libraries (mosaic-app PUBLIC mosaic)
endif()

add_library (treetops SHARED src/treetops.cpp)
target_link_libraries (treetops gdal sqlite3 spatialite raster util geotools)
add_executable (treetops-app src/apps/treetops.cpp)
target_link_libraries (treetops-app PUBLIC treetops)


add_library (spectral SHARED src/spectral.cpp)
target_link_libraries (spectral gdal sqlite3 spatialite raster util geotools)
add_executable (spectral-app src/apps/spectral.cpp)
target_link_libraries (spectral-app PUBLIC spectral)


add_library(pointstats SHARED src/pointstats.cpp src/lasutil.cpp src/pointstream.cpp)
target_link_libraries (pointstats gdal geos CGAL SFCGAL laszip las boost_system boost_filesystem raster util geotools)
if(${WITH_GUI})
	add_executable (pointstats-app src/apps/pointstats.cpp src/apps/pointstats_ui.cpp src/apps/crs_selector_ui.cpp)
	target_link_libraries (pointstats-app PUBLIC pointstats Qt5::Gui Qt5::Core Qt5::Widgets)
else()
	add_executable (pointstats-app src/apps/pointstats.cpp)
	target_link_libraries (pointstats-app PUBLIC pointstats)
endif()

add_executable (lasboundary src/lasboundary.cpp src/lasutil.cpp)
target_link_libraries (lasboundary PUBLIC las gdal las laszip geos boost_system boost_filesystem raster util geotools)


add_library(lasvalidate SHARED src/lasvalidate.cpp src/lasutil.cpp)
target_link_libraries (lasvalidate las laszip geos raster util geotools)
add_executable (lasvalidate-app src/apps/lasvalidate.cpp)
target_link_libraries (lasvalidate-app PUBLIC lasvalidate)


add_executable (pointzonalstats src/pointzonalstats.cpp src/lasutil.cpp)
target_link_libraries (pointzonalstats gdal las laszip geos boost_system boost_filesystem raster util geotools)


add_executable (lasmerge src/lasmerge.cpp src/lasutil.cpp)
target_link_libraries (lasmerge las laszip boost_system boost_filesystem util geotools)


add_executable (lasreclass src/lasreclass.cpp src/lasutil.cpp)
target_link_libraries (lasreclass las laszip boost_system boost_filesystem util geotools)


add_executable (lasclip src/lasclip.cpp src/lasutil.cpp)
target_link_libraries (lasclip gdal las laszip geos boost_system boost_filesystem util geotools)


add_executable (las2csrs src/las2csrs.cpp src/lasutil.cpp)
target_link_libraries (las2csrs gdal las laszip proj boost_system boost_filesystem util raster geotools)



install (TARGETS flood DESTINATION bin)

install (TARGETS treetops DESTINATION lib)
install (TARGETS treetops-app DESTINATION bin)
install (FILES include/treetops.hpp DESTINATION include)

install (TARGETS spectral DESTINATION lib)
install (TARGETS spectral-app DESTINATION bin)
install (FILES include/spectra.hpp DESTINATION include)

install (TARGETS rasterstats DESTINATION bin)

install (TARGETS mosaic-app DESTINATION bin)

install (TARGETS pointstats DESTINATION lib COMPONENT libraries)
install (TARGETS pointstats-app RUNTIME DESTINATION bin COMPONENT applications)
install (FILES include/pointstats.hpp DESTINATION include COMPONENT headers)

install (TARGETS lasvalidate DESTINATION lib)
install (TARGETS lasvalidate-app DESTINATION bin)
install (FILES include/lasvalidate.hpp DESTINATION include)

install (TARGETS pointzonalstats DESTINATION bin)

install (TARGETS lasmerge DESTINATION bin)

install (TARGETS lasreclass DESTINATION bin)

install (TARGETS lasclip DESTINATION bin)

install (TARGETS las2csrs DESTINATION bin)
install (DIRECTORY share/ DESTINATION share PATTERN "share/*")




