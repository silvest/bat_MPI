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

# Standalone MPI regression tests (the legacy autotools suite predates this MPI fork).
option(BAT_BUILD_PT_TESTS "Build Parallel Tempering MPI regression tests" OFF)
if(BAT_BUILD_PT_TESTS)
  enable_testing()
  add_executable(ParallelTempering_MPI_TEST test/ParallelTempering_MPI_TEST.cxx)
  target_link_libraries(ParallelTempering_MPI_TEST PRIVATE BAT::BAT MPI::MPI_CXX)
  foreach(case IN ITEMS serial partial two_temperatures factorized_warmup disabled)
    if(case STREQUAL "serial")
      set(nproc 1)
    elseif(case STREQUAL "two_temperatures")
      set(nproc 5)
    else()
      set(nproc 3)
    endif()
    add_test(NAME PT_${case} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nproc}
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:ParallelTempering_MPI_TEST> ${MPIEXEC_POSTFLAGS} ${case})
    set_tests_properties(PT_${case} PROPERTIES TIMEOUT 120 RUN_SERIAL TRUE)
  endforeach()
endif()

if(BAT_BUILD_PT_TESTS)
  add_executable(ProposalFunction_MPI_TEST test/ProposalFunction_MPI_TEST.cxx)
  target_link_libraries(ProposalFunction_MPI_TEST PRIVATE BAT::BAT MPI::MPI_CXX)
  foreach(case IN ITEMS factorized multivariate pt)
    add_test(NAME Proposal_${case} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:ProposalFunction_MPI_TEST> ${MPIEXEC_POSTFLAGS} ${case})
    set_tests_properties(Proposal_${case} PROPERTIES TIMEOUT 120 RUN_SERIAL TRUE)
  endforeach()
endif()
