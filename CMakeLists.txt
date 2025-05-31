cmake_minimum_required(VERSION 4.0.0)

project(kursovaya)

find_program(CMAKE_CXX_COMPILER NAMES g++ clang++) 
if(CMAKE_CXX_COMPILER)
   message(STATUS "Using compiler: ${CMAKE_CXX_COMPILER}")
else()
   message(WARNING "No C++ compiler found. Using default.")
endif()

include_directories(include)

add_executable(generate src/gridGeneration.cpp src/structures.cpp src/u_functions.cpp "include/common_includes.h")

add_executable(solve src/hexFEM.cpp src/input.cpp src/structures.cpp src/u_functions.cpp "include/common_includes.h" "include/input.h")