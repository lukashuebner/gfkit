# SFKIT_BENCHMARK_DEBUG enables Backward, ASan and UBSan. We have to disable malloc_count though, as it is not
# compatible with either of those three.

if (SFKIT_BENCHMARK_ENABLE_BACKWARD_CPP)
    add_executable(sfkit-bench sfkit-bench.cpp ${BACKWARD_ENABLE})
    target_link_libraries(sfkit-bench PRIVATE Backward::Backward)
else ()
    add_executable(sfkit-bench sfkit-bench.cpp)
endif ()

if (SFKIT_BENCHMARK_ENABLE_MALLOC_COUNT)
    target_compile_definitions(sfkit-bench PRIVATE "ENABLE_MALLOC_COUNT")
endif ()

target_sources(sfkit-bench PRIVATE benchmark.cpp compress.cpp print_ds_stats.cpp)

set(MALLOC_COUNT_LIBS "$<$<BOOL:${SFKIT_BENCHMARK_ENABLE_MALLOC_COUNT}>:malloc_count;dl>")
target_link_libraries(sfkit-bench PRIVATE sfkit CLI11 Catch2::Catch2 "${MALLOC_COUNT_LIBS}")
# TODO Add -march=native flag via option() and presets target_compile_options(sfkit PRIVATE "-march=native")

if (SFKIT_BENCHMARK_ENABLE_SANITIZERS)
    add_address_sanitizer(sfkit-bench)
    add_undefined_sanitizer(sfkit-bench)
endif ()
