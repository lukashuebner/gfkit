#include <fstream>
#include <sstream>
#include <string>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <catch2/catch_approx.hpp>

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

void compress(
    std::string const& trees_file,
    std::string const& forest_file,
    std::string const& bp_forest_file,
    ResultsPrinter&    results_printer
) {
    constexpr uint16_t iteration = 0;
    constexpr bool     warmup    = false;

    auto log_time = [&results_printer,
                     &trees_file](std::string const& section, std::string const& variant, Timer::duration duration) {
        results_printer.print_timer(warmup, section, variant, trees_file, duration, iteration);
    };

    auto log_mem =
        [&results_printer,
         &trees_file](std::string const& section, std::string const& variant, MemoryUsage::Report const& report) {
            results_printer.print_memory(warmup, section, variant, trees_file, report, iteration);
        };

    sfkit::tskit::TSKitTreeSequence tree_sequence(trees_file);

    Timer       timer;
    MemoryUsage memory_usage;

    // Compress the .trees to .forest file
    memory_usage.start();
    timer.start();

    sfkit::sequence::GenomicSequenceFactory dag_sequence_factory(tree_sequence);
    sfkit::dag::DAGForestCompressor         dag_forest_compressor(tree_sequence);
    sfkit::dag::DAGCompressedForest         dag_forest   = dag_forest_compressor.compress(dag_sequence_factory);
    sfkit::sequence::GenomicSequence        dag_sequence = dag_sequence_factory.move_storage();
    do_not_optimize(dag_forest);
    do_not_optimize(dag_sequence);

    log_time("compress_forest_and_sequence", "sfkit_dag", timer.stop());
    log_mem("compress_forest_and_sequence", "sfkit_dag", memory_usage.stop());

    // Save the compressed forest and sequence to a .forest file
    memory_usage.start();
    timer.start();

    sfkit::io::CompressedForestIO::save(forest_file, dag_forest, dag_sequence);

    log_time("save_forest_file", "sfkit_dag", timer.stop());
    log_mem("save_forest_file", "sfkit_dag", memory_usage.stop());

    // Compress the .trees to .bpforest file

    memory_usage.start();
    timer.start();

    using SetOfSampleSets = sfkit::samples::SetOfSampleSets<1>;
    sfkit::bp::BPForestCompressor           bp_forest_compressor(tree_sequence);
    sfkit::sequence::GenomicSequenceFactory bp_sequence_factory(tree_sequence);
    sfkit::bp::BPCompressedForest           bp_forest   = bp_forest_compressor.compress(bp_sequence_factory);
    sfkit::sequence::GenomicSequence        bp_sequence = bp_sequence_factory.move_storage();

    log_time("compress_forest_and_sequence", "sfkit_bp", timer.stop());
    log_mem("compress_forest_and_sequence", "sfkit_bp", memory_usage.stop());

    // Save the BP-compressed forest and sequence to a .bpforest file
    memory_usage.start();
    timer.start();

    sfkit::io::CompressedForestIO::save(bp_forest_file, bp_forest, bp_sequence);

    log_time("save_forest_file", "sfkit_bp", timer.stop());
    log_mem("save_forest_file", "sfkit_bp", memory_usage.stop());
}
