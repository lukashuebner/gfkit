#include <sfkit/SuccinctForest.hpp>
#include <sfkit/tskit/tskit.hpp>

#include "BenchmarkRunner.hpp"
#include "ResultsPrinter.hpp"
#include "perf.hpp"
#include "sfkit/SuccinctForest.hpp"
#include "sfkit/bp/BPForestCompressor.hpp"
#include "sfkit/dag/DAGCompressedForest.hpp"
#include "sfkit/dag/DAGForestCompressor.hpp"
#include "sfkit/io/CompressedForestIO.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"
#include "sfkit/stats/AlleleFrequencySpectrum.hpp"
#include "sfkit/tskit/tskit.hpp"
#include "timer.hpp"

void dataset_stats(
    std::string const& trees_file,
    std::string const& dag_forest_file,
    std::string const& bp_forest_file,
    ResultsPrinter&    results_printer
) {
    constexpr bool    warmup    = false;
    constexpr uint8_t iteration = 0;

    // auto log_time = [&results_printer,
    //                  &trees_file,
    //                  warmup,
    //                  iteration](std::string const& section, std::string const& variant, Timer::duration duration) {
    //     results_printer.print_timer(false, section, variant, trees_file, duration, iteration);
    // };

    // auto log_mem = [&results_printer,
    //                 iteration,
    //                 &trees_file,
    //                 warmup](std::string const& section, std::string const& variant, MemoryUsage::Report const& report) {
    //     results_printer.print_memory(warmup, section, variant, trees_file, report, iteration);
    // };

    // auto log_stat =
    //     [&results_printer,
    //      &trees_file,
    //      warmup,
    //      iteration](std::string const& variable, std::string const& variant, auto const value, std::string const& unit) {
    //         results_printer.print(warmup, "stats", variant, trees_file, variable, value, unit, iteration);
    //     };

    // Timer           timer;
    // MemoryUsage     memory_usage;
    // BenchmarkRunner bench(memory_usage, timer, log_mem, log_time);

    // --- Load the datasets ---
    sfkit::tskit::TSKitTreeSequence tree_sequence(trees_file);
    sfkit::dag::DAGCompressedForest         dag_forest;

    sfkit::sequence::GenomicSequenceFactory dag_sequence_factory(tree_sequence);
    sfkit::sequence::GenomicSequence        dag_sequence;
    sfkit::io::CompressedForestIO::load(dag_forest_file, dag_forest, dag_sequence);

    // --- Benchmark loading the BP-based tree sequence ---
    sfkit::bp::BPCompressedForest           bp_forest;
    sfkit::sequence::GenomicSequenceFactory bp_sequence_factory(tree_sequence);
    sfkit::sequence::GenomicSequence        bp_sequence;
    sfkit::io::CompressedForestIO::load(bp_forest_file, bp_forest, bp_sequence);

    // --- Number of trees ---
    auto const num_trees = tree_sequence.num_trees();
    if (num_trees != dag_forest.num_trees() || num_trees != bp_forest.num_trees()) {
        std::cerr << "Number of trees mismatch: TSKit: " << num_trees
                  << " DAG: " << dag_forest.num_trees()
                  << " BP: " << bp_forest.num_trees() << std::endl;
    }
    log_stat("num_trees", "all", num_trees, "1");

    // --- Number of samples in the trees ---
    auto const num_samples = tree_sequence.num_samples();
    if (num_samples != dag_forest.num_samples() || num_samples != bp_forest.num_samples()) {
        std::cerr << "Number of samples mismatch: TSKit: " << num_samples
                  << " DAG: " << dag_forest.num_samples()
                  << " BP: " << bp_forest.num_samples() << std::endl;
    }
    log_stat("num_samples", "all", num_samples, "1");

    // --- Number of sites in the genomic sequence ---
    auto const num_sites = tree_sequence.num_sites();
    if (num_sites != dag_sequence.num_sites() || num_sites != bp_sequence.num_sites()) {
        std::cerr << "Number of sites mismatch: TSKit: " << num_sites
                  << " DAG: " << dag_sequence.num_sites()
                  << " BP: " << bp_sequence.num_sites() << std::endl;
    }
    log_stat("num_sites", "all", num_sites, "1");

    // --- Number of mutations ---
    auto const num_mutations = tree_sequence.num_mutations();
    if (num_mutations != dag_sequence.num_mutations() || num_mutations != bp_sequence.num_mutations()) {
        std::cerr << "Number of mutations mismatch: TSKit: " << num_mutations
                  << " DAG: " << dag_sequence.num_mutations()
                  << " BP: " << bp_sequence.num_mutations() << std::endl;
    }
    log_stat("num_mutations", "all", num_mutations, "1");

    // --- Number of unique subtrees ---
    log_stat("num_unique_subtrees", "sfkit_dag", dag_forest.num_unique_subtrees(), "1");
    log_stat("num_unique_subtrees", "sfkit_bp", bp_forest.num_unique_subtrees(), "1");

    // --- Number of subtrees with mutations on them ---
    // log_stat("num_subtrees_with_mutations", "sfkit_dag", dag_forest.num_subtrees_with_mutations(), "1");
    // log_stat("num_subtrees_with_mutations", "sfkit_bp", bp_forest.num_subtrees_with_mutations(), "1");

    // --- Number of edges in the sfkit DAG and tskit tree sequence ---
    log_stat("num_edges", "sfkit_dag", dag_forest.num_edges(), "1");
    log_stat("bp_size", "sfkit_bp", bp_forest.bp_size_bit(), "Bit");
    log_stat("is_ref_size", "sfkit_bp", bp_forest.is_ref_size_bit(), "Bit");
    log_stat("is_leaf_size", "sfkit_bp", bp_forest.is_leaf_size_bit(), "Bit");
    log_stat("references_size", "sfkit_bp", bp_forest.references_size_bit(), "Bit");
    log_stat("leaves_size", "sfkit_bp", bp_forest.leaves_size_bit(), "Bit");
    log_stat("num_edges", "tskit", tree_sequence.num_edges(), "1");
}
