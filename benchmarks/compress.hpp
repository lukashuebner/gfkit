#pragma once

#include <string>

#include "ResultsPrinter.hpp"

void compress(
    std::string const& trees_file,
    std::string const& forest_file,
    std::string const& bp_forest_file,
    ResultsPrinter&    results_printer
);
