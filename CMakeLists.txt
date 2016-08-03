cmake_minimum_required(VERSION 3.1.0)
set(CMAKE_CXX_STANDARD 11)
project(dealii_heat_problem)
include_directories("include")
# deal.II, deal2lkit, pi-DoMUS
find_package(deal.II 8.4 REQUIRED HINTS $ENV{DEAL_II_DIR})
find_package(deal2lkit 1.0 REQUIRED HINTS $ENV{DEAL2LKIT_DIR})
include_directories(include ${DEAL_II_INCLUDE_DIRS} ${D2K_INCLUDE_DIRS} $ENV{PIDOMUS_DIR}/include)
# Python binding
add_subdirectory(external/pybind11)
pybind11_add_module(dealii_heat_problem source/python_wrapper.cpp)
# Link deal.II, deal2lkit, and pi-DoMUS to the Python object.
target_link_libraries(dealii_heat_problem
    $ENV{DEAL_II_DIR}/lib/libdeal_II.so.8.5.0-pre
    $ENV{DEAL2LKIT_DIR}/lib/libdeal2lkit.so
    $ENV{PIDOMUS_DIR}/lib/libpidomus-lib.so)