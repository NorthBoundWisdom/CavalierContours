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
cavc_add_test(TEST_pline)
cavc_add_test(TEST_staticspatialindex)
