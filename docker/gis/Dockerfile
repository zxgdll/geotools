# This script prepares a debian instance for use as a GIS box.
# Software installed:
# - libtiff
# - libgeotiff
# - proj.4
# - libgeos
# - postgresql
# - postgis
# - [pgpointcloud]
# - gdal
# - pdal
# - liblas
# - [maybe others]

# The default postgres user is postgres. DB-related files are owned by that user.

# Use the latest debian image.
FROM debian:latest

# Choose a directory to install to, or just use the default
ENV install_dir /usr/local

# If a non-standard install_dir is used, it'll need to be added to the lib path. This does that.
RUN echo "${install_dir}/lib" > /etc/ld.so.conf.d/gis.conf
RUN ldconfig

# Install some necessaries
RUN apt-get clean
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -y cmake wget git
RUN apt-get install -y build-essential libc6 libc6-dev libexpat1 libexpat1-dev 
RUN apt-get install -y automake libtool
RUN apt-get install -y libpython2.7-dev python-dev
RUN apt-get install -y libjson0-dev libxml2-dev libxml2-utils xsltproc docbook-xsl docbook-mathml
RUN apt-get install -y libreadline-dev bison flex
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y libcunit1-dev
RUN apt-get install chkconfig
RUN apt-get install -y libboost-dev
RUN apt-get install -y libboost-program-options1.55-dev libboost-thread1.55-dev libboost-system1.55-dev libboost-iostreams1.55-dev libboost-filesystem1.55-dev

# Make a directory to work from and go there.
RUN mkdir -p /var/buildgis
WORKDIR /var/buildgis

# Install libtiff
RUN wget http://download.osgeo.org/libtiff/tiff-4.0.6.tar.gz
RUN tar -xzf tiff-4.0.6.tar.gz
RUN rm tiff-4.0.6.tar.gz
WORKDIR tiff-4.0.6
RUN ./autogen.sh
RUN ./configure --prefix=${install_dir}
RUN make
RUN make install
RUN ldconfig
WORKDIR ..

# Install libgeotiff
RUN wget ftp://ftp.remotesensing.org/pub/geotiff/libgeotiff/libgeotiff-1.4.0.tar.gz
RUN tar -xzf libgeotiff-1.4.0.tar.gz
RUN rm libgeotiff-1.4.0.tar.gz
RUN mkdir -p libgeotiff-1.4.0/makefiles
WORKDIR libgeotiff-1.4.0/makefiles
RUN cmake -DCMAKE_INSTALL_PREFIX=${install_dir} ..
RUN make
RUN make install
WORKDIR ../..

# Install proj.4
RUN git clone https://github.com/OSGeo/proj.4.git
RUN mkdir -p proj.4/makefiles
WORKDIR proj.4/makefiles
RUN cmake -DCMAKE_INSTALL_PREFIX=${install_dir} ..
RUN make
RUN make install
WORKDIR ../..

# Install libgeos
RUN git clone https://github.com/libgeos/libgeos.git
RUN mkdir -p libgeos/makefiles
WORKDIR libgeos/makefiles
RUN cmake -DCMAKE_INSTALL_PREFIX=${install_dir} ..
RUN make
RUN make install
WORKDIR ../..

# Install PostgreSQL
RUN git clone https://github.com/postgres/postgres.git
WORKDIR postgres
RUN ./configure --prefix=${install_dir}
RUN make
RUN make install
RUN echo "${install_dir}/pgsql/lib" >> /etc/ld.so.conf.d/gis.conf
RUN ldconfig
WORKDIR ..

# Install GDAL 1.11.*
RUN git clone https://github.com/OSGeo/gdal.git
WORKDIR gdal
RUN git checkout 1.11
WORKDIR gdal
RUN ./configure --prefix=${install_dir} --with-python --with-geos=${install_dir}/bin/geos-config --with-pg=${install_dir}/bin/pg_config --with-geotiff=${install_dir}
RUN make
RUN make install
WORKDIR ../..

# Install PostGIS
RUN git clone https://github.com/postgis/postgis.git
WORKDIR postgis
RUN ./autogen.sh
RUN ./configure --with-gdalconfig=${install_dir}/bin/gdal-config --with-projdir=${install_dir} --with-pgconfig=${install_dir}/bin/pg_config --with-geosconfig=${install_dir}/bin/geos-config
RUN make
RUN make install
RUN make comments-install
WORKDIR ..

# Install PGPointCloud
# Commented out for "atomics.h" [etc] cannot be included from front end code
#RUN git clone https://github.com/pgpointcloud/pointcloud.git
#RUN mkdir -p pointcloud/makefiles
#WORKDIR pointcloud/makefiles
#RUN cmake -DCMAKE_INSTALL_PREFIX=${install_dir} ..
#RUN make
#RUN make install
#WORKDIR ../..

# Install liblas
RUN git clone https://github.com/libLAS/libLAS.git
RUN mkdir -p libLAS/makefiles
WORKDIR libLAS/makefiles
RUN cmake -DCMAKE_INSTALL_PREFIX=${install_dir} ..
RUN make
RUN make install
WORKDIR ../..


# Install PDAL
RUN git clone https://github.com/PDAL/PDAL.git
RUN mkdir -p PDAL/makefiles
WORKDIR PDAL/makefiles
RUN cmake -DCMAKE_INSTALL_PREFIX=${install_dir} ..
RUN make
RUN make install
WORKDIR ../..

# Install the GSC geoid models into proj.4's share directory.
RUN wget http://dijital.ca/files/can_grids.tar.gz
RUN tar -xzf can_grids.tar.gz
RUN rm can_grids.tar.gz
RUN cp can_grids/* ${install_dir}/share/proj

# Install vertcs override with Canadian GEOID references in gdal share.
RUN wget http://dijital.ca/files/vertcs.override.csv
RUN cp vertcs.override.csv ${install_dir}/share/gdal

# Install PostgreSQL start script. Replace the bin path in the script.
RUN groupadd postgres
RUN useradd -g postgres postgres
RUN cp postgres/contrib/start-scripts/linux /etc/init.d/postgres
RUN sed -i '/^prefix=.*/c\prefix=/usr/local' /etc/init.d/postgres
RUN chmod +x /etc/init.d/postgres
RUN chkconfig --add postgres
RUN mkdir -p /usr/local/pgsql/data
RUN chown -R postgres:postgres /usr/local/pgsql
RUN su -c "initdb -D /usr/local/pgsql/data" - postgres

