FROM ubuntu:latest

RUN apt-get clean
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -y cmake
RUN apt-get install -y wget
RUN apt-get install -y git
RUN apt-get install -y build-essential libc6 libc6-dev libexpat1 libexpat1-dev 
RUN apt-get install -y automake libtool
RUN apt-get install -y libpython2.7-dev python-dev
RUN apt-get install -y libjson0-dev libxml2-dev libxml2-utils xsltproc docbook-xsl docbook-mathml
RUN apt-get install -y libreadline-dev bison flex

RUN wget http://download.osgeo.org/libtiff/tiff-4.0.6.tar.gz
RUN tar -xzf tiff-4.0.6.tar.gz
WORKDIR tiff-4.0.6
RUN ./autogen.sh
RUN ./configure --prefix=/opt
RUN make
RUN sudo make install
RUN sudo ldconfig
WORKDIR ..

RUN wget ftp://ftp.remotesensing.org/pub/geotiff/libgeotiff/libgeotiff-1.4.0.tar.gz
RUN tar -xzf libgeotiff-1.4.0.tar.gz
RUN mkdir -p libgeotiff-1.4.0/makefiles
WORKDIR libgeotiff-1.4.0/makefiles
RUN cmake -DCMAKE_INSTALL_PREFIX=/opt ..
RUN make
RUN sudo make install
WORKDIR ../..

RUN git clone https://github.com/OSGeo/proj.4.git
RUN mkdir -p proj.4/makefiles
WORKDIR proj.4/makefiles
RUN cmake -DCMAKE_INSTALL_PREFIX=/opt ..
RUN make
RUN sudo make install
WORKDIR ../..

RUN git clone https://github.com/libgeos/libgeos.git
RUN mkdir -p libgeos/makefiles
WORKDIR libgeos/makefiles
RUN cmake -DCMAKE_INSTALL_PREFIX=/opt ..
RUN make
RUN sudo make install
WORKDIR ../..

RUN git clone https://github.com/postgres/postgres.git
WORKDIR postgres
RUN apt-get install -y zlib1g-dev
RUN ./configure --prefix=/opt
RUN make
RUN sudo make install
RUN sudo ldconfig
WORKDIR ..

RUN git clone https://github.com/OSGeo/gdal.git
WORKDIR gdal/gdal
RUN ./configure --prefix=/opt --with-python --with-geos=/opt/bin/geos-config --with-pg=/opt/bin/pg_config --with-geotiff=/opt
RUN make
RUN sudo make install
RUN sudo ldconfig
WORKDIR ../..

RUN git clone https://github.com/postgis/postgis.git
WORKDIR postgis
RUN ./autogen.sh
RUN LDFLAGS=-lstdc++ ./configure --with-pgconfig=/opt/bin/pg_config --with-geosconfig=/opt/bin/geos-config
RUN make
RUN sudo make install
RUN sudo ldconfig
RUN sudo make comments-install
WORKDIR ..

