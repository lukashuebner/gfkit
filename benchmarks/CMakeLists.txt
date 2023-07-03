# TDT_BENCHMARK_DEBUG enables Backward, ASan and UBSan. We have to disable malloc_count though, as it is not compatible
# with either of those three.

if (TDT_BENCHMARK_DEBUG)
    add_executable(sfkit sfkit.cpp ${BACKWARD_ENABLE})
    add_backward(sfkit)
else ()
    add_executable(sfkit sfkit.cpp)
    target_compile_definitions(sfkit PRIVATE "HAS_MALLOC_COUNT")
endif ()

target_link_libraries(sfkit PRIVATE tdt CLI11 Catch2::Catch2 "$<$<NOT:$<BOOL:${TDT_BENCHMARK_DEBUG}>>:malloc_count;dl>")

if (TDT_BENCHMARK_DEBUG)
    add_address_sanitizer(sfkit)
    add_undefined_sanitizer(sfkit)
endif ()
