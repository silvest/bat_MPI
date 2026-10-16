# One object build feeds the shared/static variants; object targets are internal.
function(bat_add_library name)
  cmake_parse_arguments(LIB "" "" "HEADERS;SOURCES;LINKDEF;DEPENDENCIES" ${ARGN})
  set(dict_headers)
  foreach(header IN LISTS LIB_HEADERS)
    get_filename_component(base "${header}" NAME)
    configure_file("${header}" "${PROJECT_BINARY_DIR}/include/BAT/${base}" COPYONLY)
    list(APPEND dict_headers "BAT/${base}")
  endforeach()
  include_directories(BEFORE "${PROJECT_BINARY_DIR}/include")
  ROOT_GENERATE_DICTIONARY(G__${name} ${dict_headers} MODULE ${name} LINKDEF ${LIB_LINKDEF})
  add_library(${name}_objects OBJECT ${LIB_SOURCES} G__${name}.cxx)
  set_target_properties(${name}_objects PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_include_directories(${name}_objects PRIVATE "${CMAKE_CURRENT_BINARY_DIR}"
    "${PROJECT_BINARY_DIR}/src" "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/BAT")
  target_link_libraries(${name}_objects PRIVATE ${LIB_DEPENDENCIES})
  target_compile_options(${name}_objects PRIVATE -Wall -Wno-unknown-pragmas)
  if(BAT_BUILD_SHARED_LIBS)
    add_library(${name} SHARED $<TARGET_OBJECTS:${name}_objects>)
  else()
    add_library(${name} STATIC $<TARGET_OBJECTS:${name}_objects>)
  endif()
  add_library(BAT::${name} ALIAS ${name})
  target_include_directories(${name} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  target_link_libraries(${name} PUBLIC ${LIB_DEPENDENCIES})
  set_target_properties(${name} PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 1)
  set(install_targets ${name})
  if(BAT_BUILD_STATIC_LIBS AND BAT_BUILD_SHARED_LIBS)
    add_library(${name}_static STATIC $<TARGET_OBJECTS:${name}_objects>)
    add_library(BAT::${name}_static ALIAS ${name}_static)
    target_include_directories(${name}_static PUBLIC
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    set(static_deps ${LIB_DEPENDENCIES})
    list(TRANSFORM static_deps REPLACE "^BAT::BAT$" "BAT::BAT_static")
    target_link_libraries(${name}_static PUBLIC ${static_deps})
    set_target_properties(${name}_static PROPERTIES OUTPUT_NAME ${name})
    list(APPEND install_targets ${name}_static)
  endif()
  install(TARGETS ${install_targets} EXPORT BATTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES ${LIB_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/BAT)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/lib${name}.rootmap"
    "${CMAKE_CURRENT_BINARY_DIR}/lib${name}_rdict.pcm"
    DESTINATION ${CMAKE_INSTALL_LIBDIR} OPTIONAL)
endfunction()
