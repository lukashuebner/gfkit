cmake_minimum_required(VERSION 3.12)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Project setup
project(
    tdt
    DESCRIPTION "TODO"
    LANGUAGES CXX
)

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # folder support for IDEs
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # this has to be enabled in the main CMakeLists file
    include(CTest)

    # find Doxygen
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        if (DOXYGEN_VERSION VERSION_LESS "1.9.2")
            message(
                WARNING
                    "Doxygen must be version 1.9.2 or newer. Documentation may not be displayed correctly and CI may "
                    "not pass even if checks pass locally."
            )
        endif ()
        add_custom_target(
            docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_SOURCE_DIR}/docs/Doxyfile
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            COMMENT "Generating Documentation"
            VERBATIM
        )
    else ()
        message(STATUS "Doxygen not found, not building docs")
    endif ()
endif ()

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if (EXISTS "${LOC_PATH}")
    message(
        FATAL_ERROR
            "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build "
            "subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles."
    )
endif ()

option(TDT_WARNINGS_ARE_ERRORS OFF)
option(TDT_BUILD_EXAMPLES_BENCHMARKS_AND_TESTS OFF)
option(TDT_TESTS_DISCOVER OFF)
option(TDT_EXPENSIVE_TESTS "Enable expensive tests?" OFF)

# Enable compilation with ccache. Defaults to ON if this is the main project.
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(TDT_USE_CCACHE "Globally enable ccache." ON)
else ()
    option(TDT_USE_CCACHE "Globally enable ccache." OFF)
endif ()

if (TDT_USE_CCACHE)
    include(CCache)
endif ()

add_subdirectory(extern)

add_library(tdt_base INTERFACE)

target_include_directories(tdt_base INTERFACE include)

# Set the C++ standard to C++20.
target_compile_features(tdt_base INTERFACE cxx_std_20)

list(
    APPEND
    TDT_WARNING_FLAGS
    "-Wall"
    "-Wextra"
    "-Wconversion"
    "-Wnon-virtual-dtor"
    "-Woverloaded-virtual"
    "-Wshadow"
    "-Wsign-conversion"
    "-Wundef"
    "-Wunreachable-code"
    "-Wunused"
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    list(
        APPEND
        TDT_WARNING_FLAGS
        "-Wcast-align"
        "-Wnull-dereference"
        "-Wpedantic"
        "-Wextra-semi"
        "-Wno-gnu-zero-variadic-macro-arguments"
    )
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    list(
        APPEND
        TDT_WARNING_FLAGS
        "-Wcast-align"
        "-Wnull-dereference"
        "-Wpedantic"
        "-Wnoexcept"
        "-Wsuggest-attribute=const"
        "-Wsuggest-attribute=noreturn"
        "-Wsuggest-override"
    )
endif ()

# OFF by default.
if (TDT_WARNINGS_ARE_ERRORS)
    list(APPEND TDT_WARNING_FLAGS "-Werror")
endif ()

# Link in the OpenSSL library
find_package(OpenSSL REQUIRED)
target_link_libraries(tdt_base INTERFACE OpenSSL::SSL)

# Target for user-code
add_library(tdt INTERFACE)
add_library(tdt::tdt ALIAS tdt)
target_link_libraries(tdt INTERFACE tdt_base)

# If enabled, use exceptions, otherwise use KASSERT()
option(TDT_EXCEPTION_MODE "Use exceptions to report recoverable errors." ON)
if (TDT_EXCEPTION_MODE)
    set(KASSERT_EXCEPTION_MODE 1)
    target_compile_definitions(tdt INTERFACE -DKASSERT_EXCEPTION_MODE)
    message(STATUS "Build with exceptions enabled.")
else ()
    set(KASSERT_EXCEPTION_MODE 0)
    message(STATUS "Build with exceptions disabled. Assertions are used instead.")
endif ()

# TODO Remove the communication assertions. TODO Unify the assertion levels with the C preprocessor definitions The
# assertion level controls which assertions are enabled during runtime:
#
# * Level 0: Disable all assertions
# * Level 10: Exception assertions = only enable exceptions (if not in exception mode)
# * Level 20: Light assertions = assertions that do not affect the running time of library operations significantly.
# * Level 30: Normal assertions = assertions that might slow down some operations of the library by a constant factor.
#   Should only be used in debug mode.
# * Level 60: Heavy assertions = assertions that introduce overhead which renders some library operations infeasible
#   when invoked with any significant work load.
#
# Assertion levels can be set explicitly using the -DTDT_ASSERTION_LEVEL=... flag. If no level is set explicitly, we set
# it to 10 (exceptions only) in Release mode and 30 (up to normal assertions) in Debug mode.
set(TDT_ASSERTION_LEVEL
    $<IF:$<CONFIG:Debug>,"normal","none">
    CACHE STRING "Assertion level"
)
set_property(
    CACHE TDT_ASSERTION_LEVEL
    PROPERTY STRINGS
             none
             exceptions
             light
             normal
             heavy
)
message(STATUS "Assertion level: ${TDT_ASSERTION_LEVEL}")

# If TDT_ASSERTION_LEVEL defaults to the generator expression, ${TDT_ASSERTION_LEVEL} may not be quoted However, if it
# is explicitly set to some constant string, it must be quoted Thus, all levels are listed twice, once with and without
# quotes @todo find a better solution for this problem
string(
    CONCAT KASSERT_ASSERTION_LEVEL
           $<$<STREQUAL:${TDT_ASSERTION_LEVEL},"none">:0>
           $<$<STREQUAL:"${TDT_ASSERTION_LEVEL}","none">:0>
           $<$<STREQUAL:${TDT_ASSERTION_LEVEL},"exceptions">:10>
           $<$<STREQUAL:"${TDT_ASSERTION_LEVEL}","exceptions">:10>
           $<$<STREQUAL:${TDT_ASSERTION_LEVEL},"light">:20>
           $<$<STREQUAL:"${TDT_ASSERTION_LEVEL}","light">:20>
           $<$<STREQUAL:${TDT_ASSERTION_LEVEL},"normal">:30>
           $<$<STREQUAL:"${TDT_ASSERTION_LEVEL}","normal">:30>
           $<$<STREQUAL:${TDT_ASSERTION_LEVEL},"heavy">:60>
           $<$<STREQUAL:"${TDT_ASSERTION_LEVEL}","heavy">:60>
)

# Add the KaSSERT assertion library
add_subdirectory("extern/kassert")
target_include_directories(tdt_base INTERFACE "extern/kassert/include")
target_link_libraries(tdt_base INTERFACE kassert)

# Add the tskit (succinct tree sequences) libary
add_subdirectory("extern/tskit/c" "extern/tskit")

# Add Catch2 (testing and benchmarking)
add_subdirectory("extern/Catch2")

# Add backward-cpp for pretty stack traces
add_subdirectory("extern/backward-cpp")

# Add nanobench for benchmarking
add_subdirectory("extern/nanobench")

# Add CLI11 for command line parsing
add_subdirectory("extern/CLI11")

# Add fmt for std::format for older libstdc++
add_subdirectory("extern/fmt")

# Add pasta::bit_vector for bitsets
add_subdirectory("extern/bit_vector")

# Google SparseHash for fast hash tables TODO Add as submodule

# Add xxHash for fast non-cryptographic hashing
set(BUILD_SHARED_LIBS OFF)
set(XXHASH_BUILD_ENABLE_INLINE_API ON)
set(XXHASH_BUILD_XXHSUM OFF)
add_subdirectory("extern/xxHash/cmake_unofficial" EXCLUDE_FROM_ALL)
# TODO Clean up, which libraries do we really need?
target_link_libraries(tdt_base INTERFACE xxhash pasta_bit_vector fmt)

# Tests, benchamrks, and examples are only built if this is the main project or if
# TDT_BUILD_EXAMPLES_BENCHMARKS_AND_TESTS is set (OFF by default)
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR TDT_BUILD_EXAMPLES_BENCHMARKS_AND_TESTS)
    add_subdirectory(examples)
    add_subdirectory(tests)
    add_subdirectory(benchmarks)
endif ()
