cmake_minimum_required (VERSION 3.1.0)
set (CMAKE_CXX_STANDARD 11)

project (geotools)

link_directories(/usr/local/lib)

add_executable (flood src/flood.cpp)
target_include_directories (flood PUBLIC include /usr/local/include /usr/include/gdal)
target_link_libraries (flood gdal)

add_executable (rasterstats src/rasterstats.cpp)
target_include_directories (rasterstats PUBLIC include /usr/local/include /usr/include/gdal)
target_link_libraries (rasterstats gdal)

add_executable (lasgrid src/lasgrid.cpp)
target_include_directories (lasgrid PUBLIC include /usr/local/include /usr/include/gdal)
target_link_libraries (lasgrid gdal las laszip boost_system boost_filesystem)

add_library(lasvalidate SHARED src/lasvalidate.cpp)
target_include_directories (lasvalidate PUBLIC include /usr/local/include)
target_link_libraries (lasvalidate las laszip geos)

add_executable (lasvalidate-app src/apps/lasvalidate.cpp)
target_include_directories (lasvalidate-app PUBLIC include /usr/local/include)
target_link_libraries (lasvalidate-app lasvalidate las geos)

add_executable (las2csrs src/las2csrs.cpp)
target_include_directories (las2csrs PUBLIC include /usr/include/gdal /usr/local/include)
target_link_libraries (las2csrs gdal las laszip proj boost_system boost_filesystem)

add_executable (lasstats src/lasstats.cpp)
target_include_directories (lasstats PUBLIC include /usr/local/include /usr/include/gdal)
target_link_libraries (lasstats gdal las laszip geos boost_system boost_filesystem)

add_executable (lasmerge src/lasmerge.cpp)
target_include_directories (lasmerge PUBLIC include /usr/local/include)
target_link_libraries (lasmerge las laszip boost_system boost_filesystem)

add_executable (lasreclass src/lasreclass.cpp)
target_include_directories (lasreclass PUBLIC include /usr/local/include)
target_link_libraries (lasreclass las laszip boost_system boost_filesystem)

add_executable (lasclip src/lasclip.cpp)
target_include_directories (lasclip PUBLIC include /usr/include/gdal /usr/local/include)
target_link_libraries (lasclip gdal las laszip geos boost_system boost_filesystem)

IF (OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
ENDIF ()

add_executable (mosaic src/mosaic.cpp)
target_include_directories (mosaic PUBLIC include /usr/include/gdal /usr/local/include)
target_link_libraries (mosaic gdal boost_system boost_filesystem)


install (TARGETS lasgrid DESTINATION bin)
install (TARGETS lasvalidate DESTINATION lib)
install (TARGETS lasvalidate-app DESTINATION bin)
install (FILES include/lasvalidate.hpp DESTINATION include)
install (TARGETS las2csrs DESTINATION bin)
install (DIRECTORY share/ DESTINATION share PATTERN "share/*")
install (TARGETS lasstats DESTINATION bin)
install (TARGETS lasmerge DESTINATION bin)
install (TARGETS lasreclass DESTINATION bin)
install (TARGETS lasclip DESTINATION bin)
install (TARGETS mosaic DESTINATION bin)
install (TARGETS rasterstats DESTINATION bin)

