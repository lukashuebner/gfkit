add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark PRIVATE tdt nanobench tskit cereal CLI11)
target_compile_definitions(benchmark PRIVATE)
