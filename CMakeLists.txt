cmake_minimum_required(VERSION 3.15)

project(2d_wave)

# find MPI package
find_package(MPI REQUIRED)

# Bring the headers
include_directories(inc)

# source file
file(GLOB SOURCES "src/*.cpp")

add_executable(main ${SOURCES})

target_link_libraries(main PUBLIC MPI::MPI_CXX)

