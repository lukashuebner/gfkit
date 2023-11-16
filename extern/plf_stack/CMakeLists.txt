cmake_minimum_required(VERSION 3.13)

project(
    plf_stack
    DESCRIPTION "A faster stack implementation"
    LANGUAGES C
)

add_library(plf_stack STATIC)
target_compile_features(plf_stack PRIVATE c_std_99)
set_target_properties(plf_stack PROPERTIES LINKER_LANGUAGE C)
target_include_directories(plf_stack PUBLIC ".")
target_sources(
    plf_stack
    PUBLIC plf_stack.h
)

# set(CMAKE_INCLUDE_CURRENT_DIR ON)

