#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include "ResultsPrinter.hpp"

void benchmark(
    bool const         warmup,
    uint8_t const      iteration,
    std::string const& trees_file,
    std::string const& forest_file,
    std::string const& bp_forest_file,
    ResultsPrinter&    results_printer
);
