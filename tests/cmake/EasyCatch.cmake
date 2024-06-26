include(Catch)

function (easycatch_set_kassert_flags TARGET_NAME)
    cmake_parse_arguments("SFKIT" "NO_EXCEPTION_MODE" "" "" ${ARGN})

    # Use global assertion level
    target_compile_definitions(${TARGET_NAME} PRIVATE -DKASSERT_ASSERTION_LEVEL=${KASSERT_ASSERTION_LEVEL})

    # Explicitly specify exception mode for tests, default to no exception mode
    if (NOT SFKIT_NO_EXCEPTION_MODE)
        target_compile_definitions(${TARGET_NAME} PRIVATE -DKASSERT_EXCEPTION_MODE)
    endif ()
endfunction ()

# Convenience wrapper for adding tests. This creates the target, links Catch2 and sfkit, enables warnings, and registers
# the test
#
# TARGET_NAME the target name FILES the files of the target
#
# example: register_test(mytarget FILES mytarget.cpp)
function (register_test TARGET_NAME)
    cmake_parse_arguments("SFKIT" "NO_EXCEPTION_MODE" "" "FILES;LIBRARIES" ${ARGN})
    add_executable(${TARGET_NAME} ${BACKWARD_ENABLE} ${SFKIT_FILES})
    target_link_libraries(${TARGET_NAME} PRIVATE Backward::Backward)
    target_link_libraries(${TARGET_NAME} PRIVATE Catch2::Catch2WithMain sfkit::sfkit "${SFKIT_LIBRARIES}")
    target_compile_options(${TARGET_NAME} PRIVATE ${SFKIT_WARNING_FLAGS})
    catch_discover_tests(${TARGET_NAME} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    easycatch_set_kassert_flags(${TARGET_NAME} ${ARGN})

    if (SFKIT_TESTS_ENABLE_SANITIZERS)
        add_address_sanitizer(${TARGET_NAME})
        add_undefined_sanitizer(${TARGET_NAME})
    endif ()
endfunction ()
