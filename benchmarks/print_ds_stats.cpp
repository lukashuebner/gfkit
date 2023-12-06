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
) {
    // Number of trees
    results_printer
        .print(warmup, "stats", "sfkit_dag", trees_file, "num_trees", dag_forest.num_trees(), "1", iteration);
    results_printer.print(warmup, "stats", "tskit", trees_file, "num_trees", tree_sequence.num_trees(), "1", iteration);

    // Number of unique subtrees
    results_printer.print(
        warmup,
        "stats",
        "sfkit_dag",
        trees_file,
        "num_unique_subtrees",
        dag_forest.num_unique_subtrees(),
        "1",
        iteration
    );

    // Number of samples in the trees
    results_printer
        .print(warmup, "stats", "sfkit_dag", trees_file, "num_samples", dag_forest.num_samples(), "1", iteration);
    results_printer
        .print(warmup, "stats", "tskit", trees_file, "num_samples", tree_sequence.num_samples(), "1", iteration);

    // Number of sites in the genomic sequence
    results_printer
        .print(warmup, "stats", "sfkit_dag", trees_file, "num_sites", dag_forest.num_sites(), "1", iteration);
    results_printer.print(warmup, "stats", "tskit", trees_file, "num_sites", tree_sequence.num_sites(), "1", iteration);

    // Number of mutations
    results_printer
        .print(warmup, "stats", "sfkit_dag", trees_file, "num_mutations", dag_forest.num_mutations(), "1", iteration);
    results_printer
        .print(warmup, "stats", "tskit", trees_file, "num_mutations", tree_sequence.num_mutations(), "1", iteration);

    // Number of subtrees with mutations on them
    results_printer.print(
        warmup,
        "stats",
        "sfkit_dag",
        trees_file,
        "num_subtrees_with_mutations",
        dag_forest.num_subtrees_with_mutations(),
        "1",
        iteration
    );

    // Number of edges in the DAG
    results_printer.print(warmup, "stats", "tskit", trees_file, "num_edges", tree_sequence.num_edges(), "1", iteration);

    // TODO Output CSV line
    std::cout << "# uniq subtrees | BP: " << bp_forest.num_unique_subtrees()
              << " DAG: " << dag_forest.num_unique_subtrees() << std::endl;
    // std::cout << "# backrefs | BP: " << bp_forest.num_backrefs() << std::endl;
    // std::cout << "# nodes | BP: " << bp_forest.num_nodes() << " DAG: " << dag_forest.num_nodes() << std::endl;
    // std::cout << "# leaves | BP: " << bp_forest.num_leaves() << " DAG: " << dag_forest.num_leaves() << std::endl;
    std::cout << "# samples | BP: " << bp_forest.num_samples() << " DAG: " << dag_forest.num_samples() << std::endl;
    std::cout << "# trees | BP: " << bp_forest.num_trees() << " DAG: " << dag_forest.num_trees()
              << " TS: " << tree_sequence.num_trees() << std::endl;
}
