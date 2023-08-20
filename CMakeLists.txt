cmake_minimum_required(VERSION 3.22)
project(cMathKernel)

set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")

# 等价于gcc -I
include_directories(${PROJECT_SOURCE_DIR}/include)

aux_source_directory(./src/factor SRC_FACTOR)
aux_source_directory(./src/matlab SRC_MATLAB)

add_executable(cMathKernel ./src/main.cpp ${SRC_FACTOR} ${SRC_MATLAB} src/data_gen.cpp include/data_gen.h)

#add_executable(demo demo.cpp) # 生成可执行文件
#add_library(common STATIC util.cpp) # 生成静态库
#add_library(common SHARED util.cpp) # 生成动态库或共享库