#!/bin/sh

# check if OpenMPI is cached from previous build
if [ -f "openmpi/bin/mpirun"]; then
 echo "Using cached OpenMPI"
else
 echo "Downloading OpenMPI source"
 wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.2.tar.gz
 tar xfz openmpi-4.0.2.tar.gz
 rm openmpi-4.0.2.tar.gz
 echo "Configuring and building openmpi"
 cd openmpi-4.0.2
 ./configure --prefix=`pwd`/../openmpi
 make -j 4 all
 make install
 cd 
 rm -rf openmpi-4.0.2
fi

#if [ -f mpich/lib/libmpich.so ]; then
#  echo "libmpich.so found -- nothing to build."
#else
#  echo "Downloading mpich source."
#  wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
#  tar xfz mpich-3.2.tar.gz
#  rm mpich-3.2.tar.gz
#  echo "configuring and building mpich."
#  cd mpich-3.2
#  ./configure \
#          --prefix=`pwd`/../mpich \
#          --enable-static=false \
#          --enable-alloca=true \
#          --disable-long-double \
#          --enable-threads=single \
#          --enable-fortran=no \
#          --enable-fast=all \
#          --enable-g=none \
#          --enable-timing=none
#  make -j4
#  make install
#  cd -
#  rm -rf mpich-3.2
#fi
