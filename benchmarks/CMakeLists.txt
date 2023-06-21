add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark PRIVATE tdt nanobench tskit cereal CLI11 malloc_count Catch2::Catch2)
target_compile_definitions(benchmark PRIVATE)
