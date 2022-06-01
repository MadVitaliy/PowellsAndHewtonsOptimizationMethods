cmake_minimum_required(VERSION 3.20)
project(PowellsMethod)

set(CMAKE_CXX_STANDARD 17)

add_executable(PowellsMethod main.cpp point.h point.cpp Optimization.cpp Optimization.h)
