cmake_minimum_required(VERSION 3.2)

# name of the project + language
project(2d_wave CXX)

# set cpp standard
set(CMAKE_CXX_STANDARD 11)  # C++ 11
set(CMAKE_POSITION_INDEPENDENT_CODE ON)	# initialize the POSITION_INDEPENDENT_CODE property on all target/libraries. 

# have unit tets
enable_testing()

# Bring the headers
include_directories(inc)

# source file
file(GLOB SOURCES "src/*.cpp")

# executable file
add_executable(main ${SOURCES})

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


# MPI---------------------------------------------------------------------------------------

# find MPI package
#find_package(MPI REQUIRED)
#if(MPI_NOFOUND)
#	message(ERROR "MPI is not found. Please check FindMPI.")
#endif(MPI_NOFOUND)

include_directories(/home/travis/build/ShiqiHe000/DG_wave_c/openmpi/include)
find_library(MPI NAMES openmpi HINTS /home/travis/build/ShiqiHe000/DG_wave_c/openmpi/lib)
if(MPI_NOTFOUND)
	message(ERROR "openmpi not found.")
endif(MPI_NOTFOUND)
#-----------------------------------------------------------------------------------------------

# libraires links
target_link_libraries(main PUBLIC ${MPI} )
#target_link_libraries(main PUBLIC MPI::MPI_CXX ${BLAS})

