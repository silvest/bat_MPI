cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# This CMake build currently covers the core BAT library only (BAT/ headers
# + src/). models/, test/, examples/ and doc/ are still built via autotools
# (see configure.ac / Makefile.am), and CUBA, RooStats and OpenMP support are
# not yet available through this CMake build.
project(BAT VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(BAT_BUILD_SHARED_LIBS "Build BAT as a shared library" ON)
option(BAT_BUILD_STATIC_LIBS "Build BAT as a static library" ON)

include(GNUInstallDirs)

find_package(ROOT 5.34.19 REQUIRED
  COMPONENTS Core Imt RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript
             Matrix Physics MathCore Minuit MathMore
)
include(${ROOT_USE_FILE})

find_package(MPI REQUIRED)

add_subdirectory(src)
