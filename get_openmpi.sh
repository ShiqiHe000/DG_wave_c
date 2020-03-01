#!/bin/sh

# check if OpenMPI is cached from previous build
if [ -f openmpi/lib/libmpi.so ]; then
 echo "Using cached OpenMPI"
else
 echo "Downloading OpenMPI source"
 mkdir openmpi
 wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.2.tar.gz
 tar xfz openmpi-4.0.2.tar.gz
 rm openmpi-4.0.2.tar.gz
 echo "Configuring and building openmpi"
 cd openmpi-4.0.2
 echo `pwd`
 ./configure --prefix=`pwd`/../openmpi
 make -j 4 all
 make install
 cd ..
fi

#test -n $CC && unset CC
#test -n $CXX && unset CXX
