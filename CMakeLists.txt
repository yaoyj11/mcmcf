cmake_minimum_required(VERSION 3.10)
project(mcmcf)

set(CMAKE_CXX_STANDARD 11)

add_executable(mcmcf main.cc fractional_packing.cc fractional_packing.h flow_solution.cc flow_solution.h demand.h edge.h)