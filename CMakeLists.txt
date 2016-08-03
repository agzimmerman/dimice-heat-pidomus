cmake_minimum_required(VERSION 3.1.0)
set(CMAKE_CXX_STANDARD 11)
project(dealii_heat_problem)
include_directories("include")
# deal.II, deal2lkit, pi-DoMUS
find_package(deal.II 8.4 REQUIRED HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR})
find_package(deal2lkit 1.0 REQUIRED HINTS ${DEAL2LKIT_DIR} $ENV{DEAL2LKIT_DIR})
include_directories(include ${DEAL_II_INCLUDE_DIRS} ${D2K_INCLUDE_DIRS} $ENV{PIDOMUS_DIR}/include)
# Python binding
add_subdirectory(external/pybind11)
pybind11_add_module(dealii_heat_problem source/python_wrapper.cpp)