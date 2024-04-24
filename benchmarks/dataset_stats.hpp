#pragma once

#include <sfkit/SuccinctForest.hpp>
#include <sfkit/tskit/tskit.hpp>

#include "ResultsPrinter.hpp"

void dataset_stats(
    std::string const& trees_file,
    std::string const& dag_forest_file,
    std::string const& bp_forest_file,
    ResultsPrinter&    results_printer
);
