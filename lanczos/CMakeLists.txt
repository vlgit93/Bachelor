cmake_minimum_required(VERSION 2.8)
PROJECT(lanczos)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_EXPORT_COMPILE_COMMANDS On)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -lgsl -lgslcblas -lboost_program_options -march=native")

ADD_EXECUTABLE(lanczos lanczos.cpp)


