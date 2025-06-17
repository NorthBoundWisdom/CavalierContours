function(AddCoverage target)
    if(NOT ADD_COVERAGE OR WIN32 OR APPLE)
        return()
    endif()

    # Parse optional arguments for coverage directories and files
    set(options "")
    set(oneValueArgs "")
    set(multiValueArgs COVERAGE_DIRS COVERAGE_FILES)
    cmake_parse_arguments(COV "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        add_gcc_coverage(${target} COVERAGE_DIRS ${COV_COVERAGE_DIRS} COVERAGE_FILES ${COV_COVERAGE_FILES})
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_clang_coverage(${target} COVERAGE_DIRS ${COV_COVERAGE_DIRS} COVERAGE_FILES ${COV_COVERAGE_FILES})
    else()
        message(STATUS "Coverage not supported for compiler: ${CMAKE_CXX_COMPILER_ID}")
    endif()

endfunction()


# Add compile and link options for coverage to libraries
function(add_coverage_compile_options target)
    if(NOT ADD_COVERAGE OR WIN32 OR APPLE)
        return()
    endif()

    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(STATUS "Added GCC coverage flags to ${target}")
        target_compile_options(${target} PRIVATE --coverage -O0 -fno-inline -fno-eliminate-unused-debug-types -fprofile-arcs -ftest-coverage)
        target_link_options(${target} PRIVATE --coverage -lgcov)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(STATUS "Added Clang coverage flags to ${target}")
        target_compile_options(${target} PRIVATE -fprofile-instr-generate -fcoverage-mapping -O0 -fno-inline)
        target_link_options(${target} PRIVATE -fprofile-instr-generate -fcoverage-mapping)
    endif()
endfunction()

function(add_gcc_coverage target)
    # Parse optional arguments for coverage directories and files
    set(options "")
    set(oneValueArgs "")
    set(multiValueArgs COVERAGE_DIRS COVERAGE_FILES)
    cmake_parse_arguments(COV "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    find_program(LCOV_PATH lcov)
    if(NOT LCOV_PATH)
        message(WARNING "lcov not found! Coverage will not be available.")
        return()
    endif()

    find_program(GENHTML_PATH genhtml)
    if(NOT GENHTML_PATH)
        message(WARNING "genhtml not found! Coverage will not be available.")
        return()
    endif()

    message(STATUS "GCC Coverage enabled for ${target}")

    set(target_name Coverage_${target})

    # Build filter commands based on coverage directories and files
    set(EXTRACT_COMMANDS "")
    set(COMBINE_ARGS "")

    if(COV_COVERAGE_FILES)
        # If specific files are specified, ONLY include those files
        foreach(file ${COV_COVERAGE_FILES})
            # Create a safe filename from file path
            get_filename_component(filename "${file}" NAME_WE)
            string(REPLACE "/" "_" safe_file_name "${filename}")
            list(APPEND EXTRACT_COMMANDS
                "COMMAND" "${LCOV_PATH}" "-e" "${target_name}.info" "${file}" "-o" "${safe_file_name}_temp.info" "--rc" "branch_coverage=1" "--ignore-errors" "unused")
            list(APPEND COMBINE_ARGS "${safe_file_name}_temp.info")
        endforeach()

        set(FILTER_COMMANDS
            ${EXTRACT_COMMANDS}
            COMMAND ${LCOV_PATH} -a ${COMBINE_ARGS} -o filtered.info --rc branch_coverage=1
            COMMAND rm -f ${COMBINE_ARGS}
        )
        message(STATUS "GCC Coverage will ONLY include files: ${COV_COVERAGE_FILES}")
    elseif(COV_COVERAGE_DIRS)
        # Process directories only if no specific files are specified
        foreach(dir ${COV_COVERAGE_DIRS})
            # Create a safe filename from directory path
            string(REPLACE "/" "_" safe_dir_name "${dir}")
            list(APPEND EXTRACT_COMMANDS
                "COMMAND" "${LCOV_PATH}" "-e" "${target_name}.info" "${CMAKE_SOURCE_DIR}/${dir}/*" "-o" "${safe_dir_name}_temp.info" "--rc" "branch_coverage=1" "--ignore-errors" "unused")
            list(APPEND COMBINE_ARGS "${safe_dir_name}_temp.info")
        endforeach()

        set(FILTER_COMMANDS
            ${EXTRACT_COMMANDS}
            COMMAND ${LCOV_PATH} -a ${COMBINE_ARGS} -o combined.info --rc branch_coverage=1
            COMMAND ${LCOV_PATH} -r combined.info '/usr/include/*' -o filtered.info --rc branch_coverage=1 --ignore-errors unused
            COMMAND rm -f ${COMBINE_ARGS} combined.info
        )
        message(STATUS "GCC Coverage will focus on directories: ${COV_COVERAGE_DIRS}")
    else()
        # Default behavior: exclude system headers and third-party code
        set(FILTER_COMMANDS
            COMMAND ${LCOV_PATH} -r ${target_name}.info '/usr/include/*' -o filtered.info --rc branch_coverage=1 --ignore-errors unused
        )
    endif()

    add_custom_target(${target_name}
        COMMENT "Running GCC coverage for ${target}..."
        COMMAND ${LCOV_PATH} -d . --zerocounters
        COMMAND $<TARGET_FILE:${target}>
        COMMAND ${LCOV_PATH} -d . --capture -o ${target_name}.info --ignore-errors inconsistent,usage,version,mismatch --rc branch_coverage=1
        ${FILTER_COMMANDS}
        COMMAND ${GENHTML_PATH} -o ${target_name} filtered.info --legend --ignore-errors inconsistent --branch-coverage --rc branch_coverage=1
        COMMAND rm -rf ${target_name}.info filtered.info
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endfunction()

function(add_clang_coverage target)
    # Parse optional arguments for coverage directories and files
    set(options "")
    set(oneValueArgs "")
    set(multiValueArgs COVERAGE_DIRS COVERAGE_FILES)
    cmake_parse_arguments(COV "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    find_program(LLVM_COV_PATH llvm-cov)
    if(NOT LLVM_COV_PATH)
        message(WARNING "llvm-cov not found! Coverage will not be available.")
        return()
    endif()

    find_program(LLVM_PROFDATA_PATH llvm-profdata)
    if(NOT LLVM_PROFDATA_PATH)
        message(WARNING "llvm-profdata not found! Coverage will not be available.")
        return()
    endif()

    message(STATUS "Clang Coverage enabled for ${target}")

    # Build filter arguments for llvm-cov
    set(FILTER_ARGS "")
    set(SOURCES_LIST "")

    if(COV_COVERAGE_FILES)
        # If specific files are specified, ONLY include those files
        foreach(file ${COV_COVERAGE_FILES})
            # Add the file directly (it should be the full path already)
            list(APPEND SOURCES_LIST "${file}")
        endforeach()

        if(SOURCES_LIST)
            list(APPEND FILTER_ARGS "--sources")
            list(APPEND FILTER_ARGS ${SOURCES_LIST})
        endif()
        message(STATUS "Clang Coverage will ONLY include files: ${COV_COVERAGE_FILES}")
    elseif(COV_COVERAGE_DIRS)
        # Process directories only if no specific files are specified
        foreach(dir ${COV_COVERAGE_DIRS})
            # Find all source files in the specified directory
            file(GLOB_RECURSE DIR_SOURCES "${CMAKE_SOURCE_DIR}/${dir}/*.cpp" "${CMAKE_SOURCE_DIR}/${dir}/*.c" "${CMAKE_SOURCE_DIR}/${dir}/*.h" "${CMAKE_SOURCE_DIR}/${dir}/*.hpp")
            list(APPEND SOURCES_LIST ${DIR_SOURCES})
        endforeach()

        if(SOURCES_LIST)
            list(APPEND FILTER_ARGS "--sources")
            list(APPEND FILTER_ARGS ${SOURCES_LIST})
        endif()
        message(STATUS "Clang Coverage will focus on directories: ${COV_COVERAGE_DIRS}")
    else()
        # Default behavior: exclude system headers and third-party code
        list(APPEND FILTER_ARGS "--ignore-filename-regex=/usr/include/.*")
    endif()

    set(target_name Coverage_${target})
    add_custom_target(${target_name}
        COMMENT "Running Clang coverage for ${target}..."
        COMMAND rm -f ${target}.profraw
        COMMAND LLVM_PROFILE_FILE=${target}.profraw $<TARGET_FILE:${target}>
        COMMAND ${LLVM_PROFDATA_PATH} merge -sparse ${target}.profraw -o ${target}.profdata
        # Generate enhanced HTML report using llvm-cov built-in options
        COMMAND ${LLVM_COV_PATH} show $<TARGET_FILE:${target}> -instr-profile=${target}.profdata
            -format=html
            -output-dir=${target_name}
            --show-branches=percent
            --show-line-counts
            --show-regions
            --show-instantiations
            --show-expansions
            --tab-size=4
            --coverage-watermark=85,50
            ${FILTER_ARGS}
        # Generate comprehensive text report
        COMMAND ${LLVM_COV_PATH} report $<TARGET_FILE:${target}> -instr-profile=${target}.profdata
            --show-branch-summary
            --show-region-summary
            --show-instantiation-summary
            ${FILTER_ARGS}
        COMMAND rm -f ${target}.profraw ${target}.profdata
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endfunction()
