cmake_minimum_required(VERSION 3.21)

# detect if CavalierContours is a subproject or not
# do not include examples, tests, and benchmarks in that case
if (NOT DEFINED PROJECT_NAME)
  set(NOT_SUBPROJECT ON)
endif()

project(CavalierContours VERSION 0.1.0)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
enable_language(C)

# Allow super-projects to override options
if (POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)
endif ()

# Build options
option(CAVC_HEADER_ONLY "C++ header only library of CavalierContours, if ON then no library is built" OFF)
option(CAVC_BUILD_SHARED_LIB "Build the C API CavalierContours dynamic shared library (SET OFF for static library)" ON)
option(CAVC_INSTALL_SDK "Install CavalierContours headers/libraries and CMake package files" OFF)
if (NOT_SUBPROJECT AND NOT CAVC_HEADER_ONLY)
    include(GoogleTest)
    enable_testing()
endif()

if(CAVC_HEADER_ONLY)
  set(CAVC_CPP_HEADER_ONLY_LIB ${PROJECT_NAME})
else()
  set(CAVC_CPP_HEADER_ONLY_LIB ${PROJECT_NAME}_header_only)
  set(CAVC_C_API_LIB ${PROJECT_NAME})
endif()

add_library(${CAVC_CPP_HEADER_ONLY_LIB} INTERFACE)
target_include_directories(${CAVC_CPP_HEADER_ONLY_LIB}
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(${CAVC_CPP_HEADER_ONLY_LIB} PROPERTIES
  EXPORT_NAME CavalierContoursHeaders)
add_library(CavalierContours::CavalierContoursHeaders ALIAS ${CAVC_CPP_HEADER_ONLY_LIB})

if(NOT CAVC_HEADER_ONLY)
  if(CAVC_BUILD_SHARED_LIB)
    add_library(${CAVC_C_API_LIB}
      SHARED src/cavaliercontours.cpp)
  else()
    add_library(${CAVC_C_API_LIB}
      STATIC src/cavaliercontours.cpp)
    target_compile_definitions(${CAVC_C_API_LIB}
      PUBLIC CAVC_STATIC_LIB)
  endif()
  # compiler warnings, note: not using /W4 on MSVC since it emits warning due to default of /W3
  # from cmake
  if(MSVC)
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
      target_compile_options(${CAVC_C_API_LIB} PRIVATE
        /clang:-Werror /clang:-Wall /clang:-Wextra /clang:-pedantic-errors /clang:-Wsign-conversion /clang:-Wconversion)
    endif()
    target_compile_options(${CAVC_C_API_LIB} PRIVATE
      /WX)
  else()
    target_compile_options(${CAVC_C_API_LIB} PRIVATE
      -Werror -Wall -Wextra -pedantic-errors -Wsign-conversion -Wconversion)
  endif()
  target_compile_definitions(${CAVC_C_API_LIB}
    PRIVATE CAVC_EXPORTS)
  target_include_directories(${CAVC_C_API_LIB}
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/c_api_include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  target_link_libraries(${CAVC_C_API_LIB}
    PRIVATE ${CAVC_CPP_HEADER_ONLY_LIB})
  set_target_properties(${CAVC_C_API_LIB} PROPERTIES
    EXPORT_NAME CavalierContoursCAPI)
  add_library(CavalierContours::CavalierContoursCAPI ALIAS ${CAVC_C_API_LIB})
endif()

if (NOT_SUBPROJECT AND NOT CAVC_HEADER_ONLY)
  # output all libraries and executables to the root build dir for simplicity in running tests with dynamic libs
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  add_subdirectory(examples)
  find_package(GTest REQUIRED)
  add_subdirectory(tests)
endif()

if(CAVC_INSTALL_SDK)
  set(CAVC_INSTALL_TARGETS ${CAVC_CPP_HEADER_ONLY_LIB})
  if(NOT CAVC_HEADER_ONLY)
    list(APPEND CAVC_INSTALL_TARGETS ${CAVC_C_API_LIB})
  endif()

  install(TARGETS ${CAVC_INSTALL_TARGETS}
    EXPORT CavalierContoursTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT CavalierContoursSDK
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT CavalierContoursSDK
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT CavalierContoursSDK
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/cavc/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cavc"
    COMPONENT CavalierContoursSDK
    FILES_MATCHING
    PATTERN "*.hpp"
  )

  if(NOT CAVC_HEADER_ONLY)
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/c_api_include/"
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
      COMPONENT CavalierContoursSDK
      FILES_MATCHING
      PATTERN "*.h"
    )
  endif()

  set(CAVC_CONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/CavalierContours")
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/CavalierContoursConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )

  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CavalierContoursConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/CavalierContoursConfig.cmake"
    INSTALL_DESTINATION "${CAVC_CONFIG_INSTALL_DIR}"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

  install(EXPORT CavalierContoursTargets
    FILE CavalierContoursTargets.cmake
    NAMESPACE CavalierContours::
    DESTINATION "${CAVC_CONFIG_INSTALL_DIR}"
    COMPONENT CavalierContoursSDK
  )

  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/CavalierContoursConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/CavalierContoursConfigVersion.cmake"
    DESTINATION "${CAVC_CONFIG_INSTALL_DIR}"
    COMPONENT CavalierContoursSDK
  )
endif()
