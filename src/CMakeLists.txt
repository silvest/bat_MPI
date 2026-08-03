# Headers in the order required for dictionary generation: cling (ROOT >= 6.08)
# needs headers listed least-dependent-first (kept in sync with src/Makefile.am).
set(BAT_DICT_HEADERS
  BAT/BCHistogramBase.h
  BAT/BCLog.h
  BAT/BCMath.h
  BAT/BCVariable.h

  BAT/BCDataPoint.h
  BAT/BCDataSet.h
  BAT/BCH1D.h
  BAT/BCH2D.h
  BAT/BCParameter.h

  BAT/BCAux.h
  BAT/BCPrior.h

  BAT/BCVariableSet.h

  BAT/BCConstantPrior.h
  BAT/BCGaussianPrior.h
  BAT/BCCombinedPrior.h
  BAT/BCSplitGaussianPrior.h
  BAT/BCCauchyPrior.h
  BAT/BCObservable.h
  BAT/BCTF1Prior.h
  BAT/BCTF1LogPrior.h
  BAT/BCTH1Prior.h
  BAT/BCPositiveDefinitePrior.h

  BAT/BCParameterSet.h

  BAT/BCEngineMCMC.h
  BAT/BCIntegrate.h
  BAT/BCModel.h
  BAT/BCEmptyModel.h
  BAT/BCPriorModel.h
  BAT/BCModelManager.h
)

# Same headers, without the "BAT/" prefix, for the public install location.
set(BAT_PUBLIC_HEADERS)
foreach(_hdr ${BAT_DICT_HEADERS})
  string(REGEX REPLACE "^BAT/" "" _hdr_name "${_hdr}")
  list(APPEND BAT_PUBLIC_HEADERS "${CMAKE_SOURCE_DIR}/BAT/${_hdr_name}")
endforeach()

set(BAT_SOURCES)
foreach(_hdr ${BAT_DICT_HEADERS})
  string(REGEX REPLACE "^BAT/(.*)\\.h$" "\\1.cxx" _src_name "${_hdr}")
  list(APPEND BAT_SOURCES "${_src_name}")
endforeach()

# ROOT_GENERATE_DICTIONARY() looks up "BAT/BC*.h" via the directories on the
# CMAKE_SOURCE_DIR include path, exactly like the tutorial in
# $ROOTSYS/tutorials/io/tree/dictionary/CMakeLists.txt.
include_directories(${CMAKE_SOURCE_DIR})

# Compute the ROOTVERSION macro the same way configure.ac / tools/build/root.m4 do:
# (major * 1000 + minor) * 1000 + patch
string(REPLACE "." ";" _bat_root_version_components "${ROOT_VERSION}")
list(GET _bat_root_version_components 0 _bat_root_major)
list(GET _bat_root_version_components 1 _bat_root_minor)
list(GET _bat_root_version_components 2 _bat_root_patch)
string(REGEX REPLACE "^0+([0-9])" "\\1" _bat_root_minor "${_bat_root_minor}")
string(REGEX REPLACE "^0+([0-9])" "\\1" _bat_root_patch "${_bat_root_patch}")
math(EXPR ROOTVERSION_NUM "(${_bat_root_major} * 1000 + ${_bat_root_minor}) * 1000 + ${_bat_root_patch}")

set(ROOTMATHMORE 1)

configure_file(${CMAKE_SOURCE_DIR}/cmake/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h @ONLY)

ROOT_GENERATE_DICTIONARY(G__BAT ${BAT_DICT_HEADERS}
  MODULE BAT
  LINKDEF LinkDef.h
)

if(NOT BAT_BUILD_SHARED_LIBS AND NOT BAT_BUILD_STATIC_LIBS)
  message(FATAL_ERROR "At least one of BAT_BUILD_SHARED_LIBS or BAT_BUILD_STATIC_LIBS must be enabled")
endif()

# Sources are compiled once (as PIC objects, so they can feed both the shared
# and the static library) and reused via target_link_libraries below.
add_library(BAT_objects OBJECT ${BAT_SOURCES} G__BAT.cxx)
set_target_properties(BAT_objects PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_compile_options(BAT_objects PRIVATE -Wall -Wno-unknown-pragmas)

target_include_directories(BAT_objects
  PRIVATE
    ${CMAKE_SOURCE_DIR}/BAT
    ${CMAKE_CURRENT_BINARY_DIR}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(BAT_objects
  PUBLIC
    ROOT::Core ROOT::Hist ROOT::Graf ROOT::Graf3d ROOT::Gpad ROOT::Tree
    ROOT::Rint ROOT::Postscript ROOT::Matrix ROOT::Physics ROOT::MathCore
    ROOT::Minuit ROOT::MathMore
  PRIVATE
    MPI::MPI_CXX
)

set(BAT_INSTALL_TARGETS)

if(BAT_BUILD_SHARED_LIBS)
  add_library(BAT SHARED)
  target_link_libraries(BAT PUBLIC BAT_objects)
  set_target_properties(BAT PROPERTIES VERSION 1.0.0 SOVERSION 1)
  if(NOT TARGET BAT::BAT)
    add_library(BAT::BAT ALIAS BAT)
  endif()
  list(APPEND BAT_INSTALL_TARGETS BAT)
endif()

if(BAT_BUILD_STATIC_LIBS)
  add_library(BAT_static STATIC)
  target_link_libraries(BAT_static PUBLIC BAT_objects)
  set_target_properties(BAT_static PROPERTIES OUTPUT_NAME BAT)
  if(NOT TARGET BAT::BAT)
    add_library(BAT::BAT ALIAS BAT_static)
  endif()
  list(APPEND BAT_INSTALL_TARGETS BAT_static)
endif()

install(TARGETS ${BAT_INSTALL_TARGETS}
  EXPORT BATTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES ${BAT_PUBLIC_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/BAT)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/libBAT.rootmap
  ${CMAKE_CURRENT_BINARY_DIR}/libBAT_rdict.pcm
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  OPTIONAL
)
