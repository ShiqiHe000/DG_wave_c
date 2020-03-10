cmake_minimum_required(VERSION 3.9)

# name of the project + language
project(2d_wave CXX)

# have unit tets
#enable_testing()

# Bring the headers
include_directories(inc)

# source file
file(GLOB SOURCES "src/*.cpp")

# executable file
add_executable(main ${SOURCES})

# set cmake compiler
set(CMAKE_CXX_COMPILER openmpi-4.0.2/bin/mpicxx)

# MPI---------------------------------------------------------------------------------------

# find MPI package
find_package(MPI REQUIRED)
if(MPI_NOFOUND)
	message(ERROR "MPI is not found. Please check FindMPI.")
endif(MPI_NOFOUND)

#include_directories(/openmpi/include)
#find_library(MPI NAMES openmpi HINTS /openmpi/lib)
#if(MPI_NOTFOUND)
#	message(ERROR "openmpi not found.")
#endif(MPI_NOTFOUND)
#-----------------------------------------------------------------------------------------------

# libraires links
target_link_libraries(main PUBLIC MPI::MPI_CXX )
#target_link_libraries(main PUBLIC MPI::MPI_CXX ${BLAS})

# OpenBLAS--------------------------------------------------------------------------------------
#include_directories(/opt/OpenBLAS/include)
#find_library(BLAS NAMES openblas HINTS /opt/OpenBLAS/lib)
#if(BLAS_NOTFOUND)	
#	message(ERROR "openblas is not found. Please check the openblas path in CMakeLists.txt")
#endif(BLAS_NOTFOUND)
#-----------------------------------------------------------------------------------------------

# silo
#include_directories(/home/shiqihe/local/silo/silo-4.10.2/include)
#find_library(SILO NAMES silo HINTS /home/shiqihe/local/silo/silo-4.10.2/lib)
#if(SILO_NOTFOUND)
#	message(ERROR "silo not found.")
#endif(SILO_NOTFOUND)

