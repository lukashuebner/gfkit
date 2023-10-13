# TDT_BENCHMARK_DEBUG enables Backward, ASan and UBSan. We have to disable malloc_count though, as it is not compatible
# with either of those three.

if (TDT_BENCHMARK_ENABLE_BACKWARD_CPP)
    add_executable(sfkit sfkit.cpp ${BACKWARD_ENABLE})
    add_backward(sfkit)
else ()
    add_executable(sfkit sfkit.cpp)
endif ()

if (TDT_BENCHMARK_ENABLE_MALLOC_COUNT)
    target_compile_definitions(sfkit PRIVATE "ENABLE_MALLOC_COUNT")
endif ()

set(MALLOC_COUNT_LIBS "$<$<BOOL:${TDT_BENCHMARK_ENABLE_MALLOC_COUNT}>:malloc_count;dl>")
target_link_libraries(sfkit PRIVATE tdt CLI11 Catch2::Catch2 "${MALLOC_COUNT_LIBS}")
# TODO Add -march=native flag via option() and presets
# target_compile_options(sfkit PRIVATE "-march=native")

if (TDT_BENCHMARK_ENABLE_SANITIZERS)
    add_address_sanitizer(sfkit)
    add_undefined_sanitizer(sfkit)
endif ()
