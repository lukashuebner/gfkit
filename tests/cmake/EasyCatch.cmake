include(Catch)

function (easycatch_set_kassert_flags TARGET_NAME)
    cmake_parse_arguments("TDT" "NO_EXCEPTION_MODE" "" "" ${ARGN})

    # Use global assertion level
    target_compile_definitions(${TARGET_NAME} PRIVATE -DKASSERT_ASSERTION_LEVEL=${KASSERT_ASSERTION_LEVEL})

    # Explicitly specify exception mode for tests, default to no exception mode
    if (NOT TDT_NO_EXCEPTION_MODE)
        target_compile_definitions(${TARGET_NAME} PRIVATE -DKASSERT_EXCEPTION_MODE)
    endif ()
endfunction ()

# Convenience wrapper for adding tests. This creates the target, links Catch2 and tdt, enables warnings, and registers
# the test
#
# TARGET_NAME the target name FILES the files of the target
#
# example: register_test(mytarget FILES mytarget.cpp)
function (register_test TARGET_NAME)
    cmake_parse_arguments("TDT" "NO_EXCEPTION_MODE" "" "FILES;LIBRARIES" ${ARGN})
    add_executable(${TARGET_NAME} ${BACKWARD_ENABLE} ${TDT_FILES})
    add_backward(${TARGET_NAME})
    target_link_libraries(
        ${TARGET_NAME} PRIVATE Catch2::Catch2WithMain tdt::tdt cereal "${TDT_LIBRARIES}" tsl::hopscotch_map
    )
    target_compile_options(${TARGET_NAME} PRIVATE ${TDT_WARNING_FLAGS})
    catch_discover_tests(${TARGET_NAME} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    easycatch_set_kassert_flags(${TARGET_NAME} ${ARGN})

    target_compile_definitions(${TARGET_NAME} PRIVATE -D_GLIBCXX_DEBUG)
    target_compile_definitions(${TARGET_NAME} PRIVATE -D_GLIBCXX_DEBUG_PEDANTIC)

    if (TDT_TEST_ENABLE_SANITIZERS)
        target_compile_options(${TARGET_NAME} PRIVATE -fsanitize=address)
        target_compile_options(${TARGET_NAME} PRIVATE -fsanitize=undefined)
        target_link_options(${TARGET_NAME} PRIVATE -fsanitize=address)
        target_link_options(${TARGET_NAME} PRIVATE -fsanitize=undefined)
    endif ()
endfunction ()
