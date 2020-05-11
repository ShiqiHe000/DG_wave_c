cmake_minimum_required(VERSION 3.9)

# name of the project + language
project(2d_wave CXX)

# specify the c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Bring the headers
include_directories(inc)

# source file
file(GLOB SOURCES "src/*.cpp")

# executable file
add_executable(main ${SOURCES})

# set compiler warning flag
set(GCC_COMPILER_FLAGS "-Wall -Wpedantic -O3 -pg")

# set cmake compiler
#set(CMAKE_CXX_COMPILER /home/travis/build/ShiqiHe000/DG_wave_c/openmpi-4.0.2/bin/mpicxx)

# MPI---------------------------------------------------------------------------------------

# find MPI package
find_package(MPI REQUIRED)
if(MPI_NOFOUND)
	message(ERROR "MPI is not found. Please check FindMPI.")
endif(MPI_NOFOUND)
#-----------------------------------------------------------------------------------------------

# libraires links
target_link_libraries(main PUBLIC MPI::MPI_CXX ${GCC_COMPILER_FLAGS})
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

