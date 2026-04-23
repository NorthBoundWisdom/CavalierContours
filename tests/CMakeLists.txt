if(CAVC_BUILD_TESTS OR CAVC_BUILD_BENCHMARKS)
    add_subdirectory(polylinefactory)
endif()

if(CAVC_BUILD_TESTS)
    add_subdirectory(tests)
endif()

if(CAVC_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
