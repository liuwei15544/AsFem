cmake_minimum_required(VERSION 3.8)
project(AsFem)

set(CMAKE_CXX_STANDARD 14)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
##############################################
### For debug mode                         ###
##############################################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O2 -Wall -Wextra")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

##############################################
### For release mode                       ###
##############################################
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -Wall -march=native -NDEBUG")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

#####################################################################
### In most cases, you don't need to modify the following lines!!!###
### In most cases, you don't need to modify the following lines!!!###
### In most cases, you don't need to modify the following lines!!!###
#####################################################################
set(SRC src/main.cpp)


add_executable(LIBs-opt ${INC} ${SRC})

