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


set(BAT_PUBLIC_HEADERS)
set(BAT_SOURCES PTDiagnostics.cxx)
foreach(header IN LISTS BAT_DICT_HEADERS)
  string(REGEX REPLACE "^BAT/(.*)\\.h$" "\\1" name "${header}")
  list(APPEND BAT_PUBLIC_HEADERS "${PROJECT_SOURCE_DIR}/${header}")
  list(APPEND BAT_SOURCES "${name}.cxx")
endforeach()

string(REPLACE "." ";" components "${ROOT_VERSION}")
list(GET components 0 major)
list(GET components 1 minor)
list(GET components 2 patch)
string(REGEX REPLACE "^0+([0-9])" "\\1" minor "${minor}")
string(REGEX REPLACE "^0+([0-9])" "\\1" patch "${patch}")
math(EXPR ROOTVERSION_NUM "(${major} * 1000 + ${minor}) * 1000 + ${patch}")
configure_file(${PROJECT_SOURCE_DIR}/cmake/config.h.in config.h @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/cmake/config_check.cxx.in config_check.cxx @ONLY)
list(APPEND BAT_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/config_check.cxx")
bat_add_library(BAT HEADERS ${BAT_PUBLIC_HEADERS} SOURCES ${BAT_SOURCES}
  LINKDEF "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h" DEPENDENCIES ${BAT_CORE_DEPENDENCIES})
