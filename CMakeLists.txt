# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(Charge_Grid Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(Charge_Grid_FORD_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/ford)
set(Charge_Grid_DOC_DIR    ${CMAKE_CURRENT_SOURCE_DIR}/doc)
set(Charge_Grid_SRC_DIR    ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}) # Just dump the test executable straight into the dir

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O -g -std=f2008 -Wall -Wextra -Wimplicit-interface -Wuse-without-only -fcheck=all -finit-real=snan" CACHE STRING "Flags used in debug compiles" FORCE)
#  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -std=f2008 -Wall -Wextra -Wimplicit-interface -Wuse-without-only" CACHE STRING "Flags used in release compiles" FORCE)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops -fstack-arrays -ffast-math -std=f2008 -Wall -Wextra -Wimplicit-interface -Wuse-without-only" CACHE STRING "Flags used in release compiles" FORCE)
endif()

# Add OMP
find_package(OpenMP)
if (OpenMP_Fortran_FOUND)
  set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS})
endif()

add_subdirectory(${Charge_Grid_SRC_DIR})
add_subdirectory(${Charge_Grid_FORD_DIR})
