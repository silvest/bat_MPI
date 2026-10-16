set(headers BCMTFAnalysisFacility.h BCMTFChannel.h BCMTFComparisonTool.h BCMTF.h
  BCMTFProcess.h BCMTFSystematic.h BCMTFSystematicVariation.h BCMTFTemplate.h)
set(sources)
foreach(header IN LISTS headers)
  string(REPLACE ".h" ".cxx" source "${header}")
  list(APPEND sources "${source}")
endforeach()
list(TRANSFORM headers PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
bat_add_library(BATmtf HEADERS ${headers} SOURCES ${sources}
  LINKDEF "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h" DEPENDENCIES BAT::BAT)
