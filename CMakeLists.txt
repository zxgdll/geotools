cmake_minimum_required(VERSION 3.0)

if (${APPLE})
	set (Qt5_DIR /usr/local/Cellar/qt5/5.6.1-1/lib/cmake/Qt5)
	set (CMAKE_C_COMPILER /usr/local/bin/gcc-6 CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/g++-6 CACHE STRING "C++ compiler" FORCE)
endif()

enable_language(C)
enable_language(CXX)

find_package(OpenMP)
find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets)

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

set (CPACK_ALL_INSTALL_TYPES Full User)

set (CPACK_COMPONENT_LIBRARIES_INSTALL_TYPES Full)
set (CPACK_COMPONENT_APPLICATIONS_INSTALL_TYPES Full User)
set (CPACK_COMPONENT_HEADERS_INSTALL_TYPES Full)

include (CPack)

project (geotools)

link_directories(/usr/local/opt/sqlite/lib /usr/local/lib /usr/lib)
include_directories(/usr/local/opt/sqlite/include include /usr/local/include /usr/include /usr/include/x86_64-linux-gnu/qt5 /usr/local/include/gdal /usr/local/include/geos /usr/include/gdal /usr/include/geos)


add_executable (flood src/flood.cpp src/geotools.cpp)
target_link_libraries (flood gdal)


add_executable (rasterstats src/rasterstats.cpp src/geotools.cpp)
target_link_libraries (rasterstats gdal)

add_executable (mosaic-app src/mosaic.cpp src/geotools.cpp)
target_link_libraries (mosaic-app gdal boost_system boost_filesystem)


add_library (treetops SHARED src/treetops.cpp src/util.cpp src/geotools.cpp)
target_link_libraries (treetops gdal sqlite3 spatialite)
add_executable (treetops-app src/apps/treetops.cpp)
target_link_libraries (treetops-app PUBLIC treetops gdal)


add_library (spectra SHARED src/spectra.cpp src/util.cpp src/geotools.cpp)
target_link_libraries (spectra gdal sqlite3 spatialite)
add_executable (spectra-app src/apps/spectra.cpp)
target_link_libraries (spectra-app PUBLIC spectra)


add_library(lasgrid SHARED src/lasgrid.cpp src/util.cpp src/lasutil.cpp src/geotools.cpp)
target_link_libraries (lasgrid las gdal las laszip boost_system boost_filesystem)
add_executable (lasgrid-app src/apps/lasgrid.cpp src/apps/lasgrid_ui.cpp)
target_link_libraries (lasgrid-app PUBLIC lasgrid Qt5::Gui Qt5::Core Qt5::Widgets)

add_executable (lasboundary src/lasboundary.cpp src/util.cpp src/lasutil.cpp src/geotools.cpp)
target_link_libraries (lasboundary PUBLIC las gdal las laszip geos boost_system boost_filesystem)

add_library(lasvalidate SHARED src/lasvalidate.cpp src/util.cpp src/lasutil.cpp src/geotools.cpp)
target_link_libraries (lasvalidate las laszip geos)
add_executable (lasvalidate-app src/apps/lasvalidate.cpp)
target_link_libraries (lasvalidate-app lasvalidate)

add_executable (lasstats src/lasstats.cpp src/util.cpp src/lasutil.cpp src/geotools.cpp)
target_link_libraries (lasstats gdal las laszip geos boost_system boost_filesystem)

add_executable (lasmerge src/lasmerge.cpp src/util.cpp src/lasutil.cpp src/geotools.cpp)
target_link_libraries (lasmerge las laszip boost_system boost_filesystem)

add_executable (lasreclass src/lasreclass.cpp src/util.cpp src/lasutil.cpp src/geotools.cpp)
target_link_libraries (lasreclass las laszip boost_system boost_filesystem)

add_executable (lasclip src/lasclip.cpp src/util.cpp src/lasutil.cpp src/geotools.cpp)
target_link_libraries (lasclip gdal las laszip geos boost_system boost_filesystem)

add_executable (las2csrs src/las2csrs.cpp src/util.cpp src/lasutil.cpp src/geotools.cpp)
target_link_libraries (las2csrs gdal las laszip proj boost_system boost_filesystem)



install (TARGETS flood DESTINATION bin)

install (TARGETS treetops DESTINATION lib)
install (TARGETS treetops-app DESTINATION bin)
install (FILES include/treetops.hpp DESTINATION include)

install (TARGETS spectra DESTINATION lib)
install (TARGETS spectra-app DESTINATION bin)
install (FILES include/spectra.hpp DESTINATION include)

install (TARGETS rasterstats DESTINATION bin)

install (TARGETS mosaic-app DESTINATION bin)

install (TARGETS lasgrid DESTINATION lib COMPONENT libraries)
install (TARGETS lasgrid-app RUNTIME DESTINATION bin COMPONENT applications)
install (FILES include/lasgrid.hpp DESTINATION include COMPONENT headers)

install (TARGETS lasvalidate DESTINATION lib)
install (TARGETS lasvalidate-app DESTINATION bin)
install (FILES include/lasvalidate.hpp DESTINATION include)

install (TARGETS lasstats DESTINATION bin)

install (TARGETS lasmerge DESTINATION bin)

install (TARGETS lasreclass DESTINATION bin)

install (TARGETS lasclip DESTINATION bin)

install (TARGETS las2csrs DESTINATION bin)
install (DIRECTORY share/ DESTINATION share PATTERN "share/*")


set (CPACK_COMPONENTS_ALL applications libraries headers)

