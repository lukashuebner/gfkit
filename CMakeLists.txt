cmake_minimum_required(VERSION 3.12)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(CCache)
include(Sanitizers)

# Project setup
project(
    sfkit
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

option(SFKIT_WARNINGS_ARE_ERRORS OFF)
option(SFKIT_BUILD_EXAMPLES_BENCHMARKS_AND_TESTS OFF) # They are also built if this is the main project.
option(SFKIT_TESTS_DISCOVER OFF)
option(SFKIT_TESTS_ENABLE_SANITIZERS "Enable sanitizers for tests?" OFF)
option(SFKIT_BENCHMARK_ENABLE_BACKWARD_CPP "Enable backward-cpp for the benchmarks" OFF)
option(SFKIT_BENCHMARK_ENABLE_MALLOC_COUNT "Enable malloc_count for the benchmarks" OFF)
option(SFKIT_BENCHMARK_ENABLE_SANITIZERS "Enable ASan and UBSan for the benchmarks" OFF)

# Enable compilation with ccache. Defaults to ON if this is the main project.
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(SFKIT_USE_CCACHE "Globally enable ccache." ON)
else ()
    option(SFKIT_USE_CCACHE "Globally enable ccache." OFF)
endif ()

if (SFKIT_USE_CCACHE)
    enable_ccache()
endif ()

add_subdirectory(extern)

add_library(sfkit_base INTERFACE)

target_include_directories(sfkit_base INTERFACE include)

# Set the C++ standard to C++20.
target_compile_features(sfkit_base INTERFACE cxx_std_20)

list(
    APPEND
    SFKIT_WARNING_FLAGS
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
        SFKIT_WARNING_FLAGS
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
        SFKIT_WARNING_FLAGS
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
if (SFKIT_WARNINGS_ARE_ERRORS)
    list(APPEND SFKIT_WARNING_FLAGS "-Werror")
endif ()

# Target for user-code
add_library(sfkit INTERFACE)
add_library(sfkit::sfkit ALIAS sfkit)
target_link_libraries(sfkit INTERFACE sfkit_base)

# If enabled, use exceptions, otherwise use KASSERT()
option(SFKIT_EXCEPTION_MODE "Use exceptions to report recoverable errors." ON)

if (SFKIT_EXCEPTION_MODE)
    set(KASSERT_EXCEPTION_MODE 1)
    target_compile_definitions(sfkit INTERFACE -DKASSERT_EXCEPTION_MODE)
    message(STATUS "Build with exceptions enabled.")
else ()
    set(KASSERT_EXCEPTION_MODE 0)
    message(STATUS "Build with exceptions disabled.  Assertions are used instead.")
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
# Assertion levels can be set explicitly using the -DSFKIT_ASSERTION_LEVEL=... flag. If no level is set explicitly, we
# set it to 10 (exceptions only) in Release mode and 30 (up to normal assertions) in Debug mode.
set(SFKIT_ASSERTION_LEVEL
    $<IF:$<CONFIG:Debug>,"normal","none">
    CACHE STRING "Assertion level"
)
set_property(CACHE SFKIT_ASSERTION_LEVEL PROPERTY STRINGS none exceptions light normal heavy)
message(STATUS "Assertion level: ${SFKIT_ASSERTION_LEVEL}")

# If SFKIT_ASSERTION_LEVEL defaults to the generator expression, ${SFKIT_ASSERTION_LEVEL} may not be quoted However, if
# it is explicitly set to some constant string, it must be quoted Thus, all levels are listed twice, once with and
# without quotes @todo find a better solution for this problem
string(
    CONCAT KASSERT_ASSERTION_LEVEL
           $<$<STREQUAL:${SFKIT_ASSERTION_LEVEL},"none">:0>
           $<$<STREQUAL:"${SFKIT_ASSERTION_LEVEL}","none">:0>
           $<$<STREQUAL:${SFKIT_ASSERTION_LEVEL},"exceptions">:10>
           $<$<STREQUAL:"${SFKIT_ASSERTION_LEVEL}","exceptions">:10>
           $<$<STREQUAL:${SFKIT_ASSERTION_LEVEL},"light">:20>
           $<$<STREQUAL:"${SFKIT_ASSERTION_LEVEL}","light">:20>
           $<$<STREQUAL:${SFKIT_ASSERTION_LEVEL},"normal">:30>
           $<$<STREQUAL:"${SFKIT_ASSERTION_LEVEL}","normal">:30>
           $<$<STREQUAL:${SFKIT_ASSERTION_LEVEL},"heavy">:60>
           $<$<STREQUAL:"${SFKIT_ASSERTION_LEVEL}","heavy">:60>
)

# Add the KaSSERT assertion library
add_subdirectory("extern/kassert")
target_include_directories(sfkit_base INTERFACE "extern/kassert/include")
target_link_libraries(sfkit_base INTERFACE kassert)

# Add the tskit (succinct tree sequences) library
add_subdirectory("extern/tskit/c" "extern/tskit")
target_link_libraries(sfkit_base INTERFACE tskit)

# Hopscotch Map as HashTable
add_subdirectory("extern/hopscotch-map")
target_link_libraries(sfkit_base INTERFACE tsl::hopscotch_map)

# Add xxHash for fast non-cryptographic hashing
set(BUILD_SHARED_LIBS OFF)
set(XXHASH_BUILD_ENABLE_INLINE_API ON)
set(XXHASH_BUILD_XXHSUM OFF)
add_subdirectory("extern/xxHash/cmake_unofficial" EXCLUDE_FROM_ALL)
target_link_libraries(sfkit_base INTERFACE xxhash)

# Add the succinct data structures library (SDSL) for bitmaps, balanced-paranthesis trees etc.
add_subdirectory("extern/sdsl")
target_link_libraries(sfkit_base INTERFACE sdsl)

# Add plf_stack for a fast stack implementation.
add_subdirectory("extern/plf_stack")
target_link_libraries(sfkit_base INTERFACE plf_stack)

# Add cereal for serialization
option(BUILD_DOC "" OFF)
option(BUILD_SANDBOX "" OFF)
option(SKIP_PERFORMANCE_COMPARISON "" ON)
add_subdirectory("extern/cereal")
target_link_libraries(sfkit_base INTERFACE cereal)

# Add Catch2 (testing and benchmarking)
add_subdirectory("extern/Catch2")

# Add backward-cpp for pretty stack traces
add_subdirectory("extern/backward-cpp")

# Add fmt for std::format for older libstdc++ TODO Remove this at least for the Release mode
add_subdirectory("extern/fmt")
target_link_libraries(sfkit_base INTERFACE fmt)

# Add CLI11 for command line parsing
add_subdirectory("extern/CLI11")

# malloc_count for analyzing memory usage
if (SFKIT_BENCHMARK_ENABLE_MALLOC_COUNT)
    add_subdirectory("extern/malloc_count")
endif ()

# Tests, benchamrks, and examples are only built if this is the main project or if
# SFKIT_BUILD_EXAMPLES_BENCHMARKS_AND_TESTS is set (OFF by default)
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR SFKIT_BUILD_EXAMPLES_BENCHMARKS_AND_TESTS)
    add_subdirectory(examples)
    add_subdirectory(tests)
    add_subdirectory(benchmarks)
endif ()
