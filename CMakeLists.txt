cmake_minimum_required(VERSION 3.22)
project(matrix_v2)

set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")

add_executable(matrix_v2 main.cpp linpack.cpp utils.cpp utils.h blas.cpp blas.h linpack.h interface.cpp interface.h)
