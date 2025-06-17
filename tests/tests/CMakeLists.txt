cmake_minimum_required(VERSION 3.21)


macro(cavc_add_test name)
    add_executable(${name} ${name}.cpp)

    target_include_directories(${name}
    PRIVATE
        ${GTEST_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/include)

    target_link_libraries(${name}
    PRIVATE
        PolylineFactory
        ${CAVC_CPP_HEADER_ONLY_LIB}
        ${CAVC_C_API_LIB}
        GTest::gtest_main
        GTest::gtest
        GTest::gmock)

    gtest_add_tests(TARGET ${name} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endmacro()

# C API tests
cavc_add_test(TEST_capi_pline)
cavc_add_test(TEST_capi_pline_function)
cavc_add_test(TEST_capi_parallel_offset)
cavc_add_test(TEST_capi_combine_plines)

# C++ API tests

cavc_add_test(TEST_plinesegment)
cavc_add_test(TEST_intrcircle2circle2)
cavc_add_test(TEST_intrlineseg2circle2)
cavc_add_test(TEST_intrlineseg2lineseg2)
cavc_add_test(TEST_pline)
cavc_add_test(TEST_staticspatialindex)

option(ADD_COVERAGE "Add coverage" 1)
include(Coverage.cmake)
if(ADD_COVERAGE)
    AddCoverage(TEST_plinesegment COVERAGE_FILES ${CMAKE_SOURCE_DIR}/include/cavc/plinesegment.hpp)
    add_coverage_compile_options(TEST_plinesegment)

    AddCoverage(TEST_intrcircle2circle2 COVERAGE_FILES ${CMAKE_SOURCE_DIR}/include/cavc/intrcircle2circle2.hpp)
    add_coverage_compile_options(TEST_intrcircle2circle2)

    AddCoverage(TEST_intrlineseg2circle2 COVERAGE_FILES ${CMAKE_SOURCE_DIR}/include/cavc/intrlineseg2circle2.hpp)
    add_coverage_compile_options(TEST_intrlineseg2circle2)

    AddCoverage(TEST_intrlineseg2lineseg2 COVERAGE_FILES ${CMAKE_SOURCE_DIR}/include/cavc/intrlineseg2lineseg2.hpp)
    add_coverage_compile_options(TEST_intrlineseg2lineseg2)

    AddCoverage(TEST_pline COVERAGE_FILES ${CMAKE_SOURCE_DIR}/include/cavc/polyline.hpp)
    add_coverage_compile_options(TEST_pline)

    AddCoverage(TEST_staticspatialindex COVERAGE_FILES ${CMAKE_SOURCE_DIR}/include/cavc/staticspatialindex.hpp)
    add_coverage_compile_options(TEST_staticspatialindex)
endif()
