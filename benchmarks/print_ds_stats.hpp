#pragma once

#include <sfkit/SuccinctForest.hpp>
#include <sfkit/tskit/tskit.hpp>

#include "ResultsPrinter.hpp"

void print_ds_stats(
    bool                                   warmup,
    std::string const&                     trees_file,
    uint16_t                               iteration,
    ResultsPrinter&                        results_printer,
    sfkit::DAGSuccinctForest const&        dag_forest,
    sfkit::BPSuccinctForest const&         bp_forest,
    sfkit::tskit::TSKitTreeSequence const& tree_sequence
);
