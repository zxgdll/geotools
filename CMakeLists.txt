cmake_minimum_required(VERSION 3.1)

set (CMAKE_POSITION_INDEPENDENT_CODE ON)

set (WITH_GUI true CACHE BOOL "If true, will build the Qt GUI.")

if(${WITH_GUI})
	add_definitions(-DWITH_GUI)
endif()

if (${APPLE})
	if(${WITH_GUI})
		set (CMAKE_PREFIX_PATH /usr/local/Cellar/qt5/5.6.1-1/lib/cmake/Qt5Widgets/)
	endif()
	set (CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE STRING "C++ compiler" FORCE)
	set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
elseif(${WIN32})
	# Qt was installed using the installer.
	# CGAL was installed from the installer, then built with cmake against Qt and OSGeo
	# Builds with VS2012 and OSGeo libs
	if(${WITH_GUI})
		set (CMAKE_PREFIX_PATH C:/Qt/5.7/msvc2013_64)
	endif()
	set (Boost_DEBUG OFF)
	set (BOOST_INCLUDEDIR C:/OSGeo4W64/include/boost-1_56)
	set (BOOST_LIBRARYDIR C:/OSGeo4W64/lib)
else()
	if(${WITH_GUI})
		set (CMAKE_PREFIX_PATH /usr/lib/x86_64-linux-gnu/cmake/Qt5)
	endif()
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -frounding-math  -D_GLIBCXX_PARALLEL")
	option(CGAL_DISABLE_ROUNDING_MATH_CHECK "Disable rounding math check in CGAL. This permits Valgrind to run." ON)
endif()

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
add_compile_options(-std=c++1y)

find_package (OpenMP)
find_package (Boost)
find_package (GDAL)
if(${WITH_GUI})
	find_package(Qt5Widgets)
endif()

set (CMAKE_AUTOUIC ON)
set (CMAKE_AUTOMOC ON)
set (CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_AUTOMOC_OPTIONS "-Iinclude")

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

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

if(${WIN32})
	link_directories(C:/Qt/5.7/msvc2013_64 C:/CGAL-4.9/makefiles/lib C:/OSGeo4W64/lib ${Boost_LIBRARY_DIR})
	include_directories(include C:/Qt/5.7/msvc2013_64/include C:/CGAL-4.9/makefiles/include C:/OSGeo4W64/include ${Boost_INCLUDE_DIR})
else()
	link_directories(${Boost_LIBRARY_DIR} /usr/local/opt/sqlite/lib /usr/local/lib /usr/lib)
	include_directories(include ${Boost_INCLUDE_DIR} ${GDAL_INCLUDE_DIR} ${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS})
endif()

add_library (geotools STATIC src/geotools.cpp)
add_library (util STATIC src/util.cpp)
if(${WIN32})
	target_link_libraries (util boost_system-vc110-mt-1_56 boost_filesystem-vc110-mt-1_56)
else()
	target_link_libraries (util boost_system boost_filesystem)
endif()
add_library (raster STATIC src/raster.cpp)
target_link_libraries (raster gdal util geotools)


add_executable (flood src/flood.cpp)
target_link_libraries (flood gdal raster util geotools)


add_executable (rasterstats src/rasterstats.cpp)
target_link_libraries (rasterstats gdal raster util geotools)


add_library (mosaic SHARED src/mosaic.cpp)
if(${WIN32})
	target_link_libraries (mosaic boost_system-vc110-mt-1_56 boost_filesystem-vc110-mt-1_56 gdal_i raster util geotools)
else()
	target_link_libraries (mosaic gdal boost_system boost_filesystem raster util geotools)
endif()

if(${WITH_GUI})
	add_executable (mosaic-app src/apps/mosaic.cpp src/apps/mosaic_ui.cpp)
	target_link_libraries (mosaic-app PUBLIC mosaic Qt5::Gui Qt5::Core Qt5::Widgets)
else()
	add_executable (mosaic-app src/apps/mosaic.cpp)
	target_link_libraries (mosaic-app PUBLIC mosaic)
endif()


add_library (treetops SHARED src/treetops.cpp)
if(${WIN32})
	target_link_libraries (treetops boost_system-vc110-mt-1_56 boost_filesystem-vc110-mt-1_56 gdal_i geos_c sqlite3_i spatialite_i raster util geotools)
else()
	target_link_libraries (treetops boost_system boost_filesystem gdal geos_c sqlite3 spatialite raster util geotools)
endif()

if(${WITH_GUI})
	add_executable (treetops-app src/apps/treetops.cpp src/apps/treetops_ui.cpp src/apps/crs_selector_ui.cpp)
	target_link_libraries (treetops-app PUBLIC treetops Qt5::Gui Qt5::Core Qt5::Widgets)
else()
	add_executable (treetops-app src/apps/treetops.cpp src/apps/treetops_ui.hpp)
	target_link_libraries (treetops-app PUBLIC treetops)
endif()


add_library (spectral SHARED src/spectral.cpp)
target_link_libraries (spectral gdal sqlite3 spatialite raster util geotools)
add_executable (spectral-app src/apps/spectral.cpp)
target_link_libraries (spectral-app PUBLIC spectral)


add_library(pointstats SHARED src/pointstats.cpp src/lasutil.cpp src/finalizedpointstream.cpp src/laspoint.cpp)
target_link_libraries (pointstats gdal geos gmp CGAL laszip las boost_system boost_filesystem raster util geotools)
if(${WITH_GUI})
	add_executable (pointstats-app src/apps/pointstats.cpp src/apps/pointstats_ui.cpp src/apps/crs_selector_ui.cpp)
	target_link_libraries (pointstats-app PUBLIC pointstats Qt5::Gui Qt5::Core Qt5::Widgets)
else()
	add_executable (pointstats-app src/apps/pointstats.cpp)
	target_link_libraries (pointstats-app PUBLIC pointstats)
endif()

add_library(pointnormalize SHARED src/pointnormalize.cpp src/pointstream.cpp src/lasutil.cpp src/laspoint.cpp)
target_link_libraries (pointnormalize gdal CGAL gmp laszip las boost_system boost_filesystem raster geotools)
if(${WITH_GUI})
	add_executable (pointnormalize-app src/apps/pointnormalize.cpp src/apps/pointnormalize_ui.cpp)
	target_link_libraries (pointnormalize-app PUBLIC pointnormalize Qt5::Gui Qt5::Core Qt5::Widgets)
else()
	add_executable (pointnormalize-app src/apps/pointnormalize.cpp)
	target_link_libraries (pointnormalize-app PUBLIC pointnormalize)
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




