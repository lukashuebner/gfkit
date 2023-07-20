#include <fstream>
#include <sstream>
#include <string>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <catch2/catch_approx.hpp>
#include <tdt/load/forest-compressor.hpp>
#include <tdt/sequence-forest.hpp>
#include <tdt/tskit.hpp>

#include "perf.hpp"
#include "tdt/graph/compressed-forest.hpp"
#include "tdt/load/compressed-forest-serialization.hpp"
#include "tdt/sequence/allele-frequency-spectrum.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"
#include "timer.hpp"

constexpr double FLOAT_EQ_EPS = 1e-4;

void compress(std::string const& trees_file, std::string const& forest_file, ResultsPrinter& results_printer) {
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

    TSKitTreeSequence tree_sequence(trees_file);

    Timer       timer;
    MemoryUsage memory_usage;

    // Compress the .trees to .forest file
    memory_usage.start();
    timer.start();

    GenomicSequenceFactory sequence_factory(tree_sequence);
    ForestCompressor       forest_compressor(tree_sequence);
    CompressedForest       forest   = forest_compressor.compress(sequence_factory);
    GenomicSequence        sequence = sequence_factory.move_storage();
    do_not_optimize(forest);
    do_not_optimize(sequence);

    log_time("compress_forest_and_sequence", "sfkit", timer.stop());
    log_mem("compress_forest_and_sequence", "sfkit", memory_usage.stop());

    // Save the compressed forest and sequence to a .forest file
    memory_usage.start();
    timer.start();

    CompressedForestIO::save(forest_file, forest, sequence);

    log_time("save_forest_file", "sfkit", timer.stop());
    log_mem("save_forest_file", "sfkit", memory_usage.stop());
}

void benchmark(
    bool const         warmup,
    uint8_t const      iteration,
    std::string const& trees_file,
    std::string const& forest_file,
    ResultsPrinter&    results_printer
) {
    auto log_time = [&results_printer, iteration, &trees_file](
                        bool const         warmup,
                        std::string const& section,
                        std::string const& variant,
                        Timer::duration    duration
                    ) {
        results_printer.print_timer(warmup, section, variant, trees_file, duration, iteration);
    };

    auto log_mem = [&results_printer, iteration, &trees_file](
                       bool const                 warmup,
                       std::string const&         section,
                       std::string const&         variant,
                       MemoryUsage::Report const& report
                   ) {
        results_printer.print_memory(warmup, section, variant, trees_file, report, iteration);
    };

    // Benchmark tree sequence loading
    Timer             timer;
    MemoryUsage       memory_usage;
    TSKitTreeSequence tree_sequence(trees_file);
    do_not_optimize(tree_sequence);
    if (!warmup) {
        log_time(warmup, "load_trees_file", "tskit", timer.stop());
    }

    // Benchmark building the DAG from the tree sequence
    memory_usage.start();
    timer.start();

    CompressedForest       forest;
    GenomicSequenceFactory sequence_factory(tree_sequence);
    GenomicSequence        sequence;

    // TODO Write IO functionality for the SequenceForest class
    CompressedForestIO::load(forest_file, forest, sequence);

    log_time(warmup, "load_forest_file", "sfkit", timer.stop());
    log_mem(warmup, "load_forest_file", "sfkit", memory_usage.stop());

    // Benchmark computing subtree sizes
    memory_usage.start();
    timer.start();

    EdgeListGraph const& dag               = forest.postorder_edges();
    auto                 num_samples_below = forest.compute_num_samples_below(forest.all_samples());
    do_not_optimize(num_samples_below);

    log_time(warmup, "compute_subtree_sizes", "sfkit", timer.stop());
    log_mem(warmup, "compute_subtree_sizes", "sfkit", memory_usage.stop());

    // Benchmark computing the AFS
    memory_usage.start();
    timer.start();

    AlleleFrequencySpectrum<PerfectDNAHasher> afs(sequence, forest);
    do_not_optimize(afs);

    log_time(warmup, "afs", "sfkit", timer.stop());
    log_mem(warmup, "afs", "sfkit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    auto reference_afs = tree_sequence.allele_frequency_spectrum();
    do_not_optimize(reference_afs);

    log_time(warmup, "afs", "tskit", timer.stop());
    log_mem(warmup, "afs", "tskit", memory_usage.stop());

    // Does our AFS match the one computed by tskit?
    bool const equal = std::ranges::equal(
        std::span(afs).subspan(1, afs.num_samples() - 1),
        std::span(reference_afs).subspan(1, afs.num_samples() - 1)
    );
    if (!equal) {
        std::cerr << "ERROR !! AFS mismatch between tskit and sf" << std::endl;
        std::exit(1);
    }

    // Benchmark computing the divergence
    // TODO The SequenceForest does not need to hold the tree_sequence at all times; it's only used during construction.
    // (Check this!)
    SequenceForest sequence_forest(std::move(tree_sequence), std::move(forest), std::move(sequence));
    SampleSet      sample_set_1(sequence_forest.num_samples());
    SampleSet      sample_set_2(sequence_forest.num_samples());
    bool           flip = false;
    for (SampleId sample: sequence_forest.all_samples()) {
        if (flip) {
            sample_set_1.add(sample);
        } else {
            sample_set_2.add(sample);
        }
        flip = !flip;
    }
    memory_usage.start();
    timer.start();

    double const sfkit_divergence = sequence_forest.divergence(sample_set_1, sample_set_2);
    do_not_optimize(sfkit_divergence);

    log_time(warmup, "divergence", "sfkit", timer.stop());
    log_mem(warmup, "divergence", "sfkit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_divergence = sequence_forest.tree_sequence().divergence(sample_set_1, sample_set_2);
    do_not_optimize(tskit_divergence);

    log_time(warmup, "divergence", "tskit", timer.stop());
    log_mem(warmup, "divergence", "tskit", memory_usage.stop());

    // Does our divergence match the one computed by tskit?
    if (sfkit_divergence != Catch::Approx(tskit_divergence).epsilon(FLOAT_EQ_EPS)) {
        std::cerr << "ERROR !! Divergence mismatch between tskit and sfkit" << std::endl;
        std::cerr << "    " << tskit_divergence << " vs. " << sfkit_divergence << " (tskit vs. sfkit)" << std::endl;
        std::exit(1);
    }

    // Benchmark computing the diversity
    memory_usage.start();
    timer.start();

    double const sfkit_diversity = sequence_forest.diversity();
    do_not_optimize(sfkit_diversity);

    log_time(warmup, "diversity", "sfkit", timer.stop());
    log_mem(warmup, "diversity", "sfkit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_diversity = sequence_forest.tree_sequence().diversity();
    do_not_optimize(tskit_diversity);

    log_time(warmup, "diversity", "tskit", timer.stop());
    log_mem(warmup, "diversity", "tskit", memory_usage.stop());

    // Does our diversity match the one computed by tskit?
    if (sfkit_diversity != Catch::Approx(tskit_diversity).epsilon(FLOAT_EQ_EPS)) {
        std::cerr << "ERROR !! Diversity mismatch between tskit and sfkit" << std::endl;
        std::cerr << "    " << tskit_diversity << " vs. " << sfkit_diversity << " (tskit vs. sfkit)" << std::endl;
        std::exit(1);
    }

    // Benchmark computing the number of segregating sites
    memory_usage.start();
    timer.start();

    uint64_t const sfkit_num_seg_sites = sequence_forest.num_segregating_sites();
    do_not_optimize(sfkit_num_seg_sites);

    log_time(warmup, "num_segregating_sites", "sfkit", timer.stop());
    log_mem(warmup, "num_segregating_sites", "sfkit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_num_seg_sites = sequence_forest.tree_sequence().diversity();
    do_not_optimize(tskit_num_seg_sites);

    log_time(warmup, "num_segregating_sites", "tskit", timer.stop());
    log_mem(warmup, "num_segregating_sites", "tskit", memory_usage.stop());

    // Does our number of segregating sites match the one computed by tskit?
    if (sfkit_num_seg_sites == Catch::Approx(tskit_num_seg_sites).epsilon(FLOAT_EQ_EPS)) {
        std::cerr << "ERROR !! Number of segregating sites mismatch between tskit and sfkit" << std::endl;
        std::cerr << "    " << tskit_num_seg_sites << " vs. " << sfkit_num_seg_sites << " (tskit vs. sfkit)"
                  << std::endl;
        std::exit(1);
    }

    // Benchmark computing Tajima's D. tskit does implement this only in Python, not in C. We're using
    // experiments/scripts/benchmark-tskits-tajimas-d.py to measure tskit's performance.
    memory_usage.start();
    timer.start();

    double const tajimas_d = sequence_forest.tajimas_d();
    do_not_optimize(tajimas_d);

    log_time(warmup, "tajimas_d", "sfkit", timer.stop());
    log_mem(warmup, "tajimas_d", "sfkit", memory_usage.stop());
    results_printer.print(warmup, "tajimas_d", "sfkit", trees_file, "tajimas_d", tajimas_d, "1", iteration);

    memory_usage.start();
    timer.start();

    // Benchmark computing the FST. tskit's implementation is in Python; see Tajima's D above.
    memory_usage.start();
    timer.start();

    double const sfkit_fst = sequence_forest.fst(sample_set_1, sample_set_2);
    do_not_optimize(sfkit_fst);

    log_time(warmup, "fst", "sfkit", timer.stop());
    log_mem(warmup, "fst", "sfkit", memory_usage.stop());
    results_printer.print(warmup, "fst", "sfkit", trees_file, "fst", sfkit_fst, "1", iteration);
}

int main(int argc, char** argv) {
    // Set-up results printer
    std::string revision              = "";
    std::string machine_id            = "";
    auto        setup_results_printer = [&revision, &machine_id]() -> ResultsPrinter {
        ResultsPrinter results_printer(std::cout, revision, machine_id);
        results_printer.print_header();
        return results_printer;
    };

    // CLI
    CLI::App app{"Forest Compression Benchmark"};

    // Compresss subcommand
    CLI::App* compress_sub =
        app.add_subcommand("compress", "Compress the tree sequence (.trees) to a succinct forest (.forest)");

    std::string trees_file = "";
    compress_sub->add_option("-f,--trees-file", trees_file, "The tree sequence input file (.trees)")
        ->check(CLI::ExistingFile)
        ->required();

    std::string forest_file = "";
    compress_sub->add_option("-i,--forest-file", forest_file, "Input file for the compressed forest")
        ->check(CLI::NonexistentPath);

    compress_sub->add_option("-r,--revision", revision, "Revision of this software (unique id, e.g. git commit hash)")
        ->default_val("undefined");

    compress_sub->add_option("-m,--machine", machine_id, "Identifier of this computer (e.g. hostname)")
        ->default_val("undefined");

    compress_sub->callback([&trees_file, &forest_file, &setup_results_printer]() {
        auto results_printer = setup_results_printer();
        std::cerr << "Compressing tree sequence " << trees_file << " -> " << forest_file << std::endl;
        compress(trees_file, forest_file, results_printer);
    });

    // Benchmark subcommand
    auto benchmark_sub = app.add_subcommand("benchmark", "Benchmark operations of the .forest and .trees files");

    // std::string trees_file = "";
    benchmark_sub
        ->add_option("-f,--trees-file", trees_file, "The tree sequence file to benchmark on. If not specified.")
        ->check(CLI::ExistingFile)
        ->required();

    // std::string forest_file = "";
    benchmark_sub->add_option("-i,--forest-file", forest_file, "Input file for the compressed forest")
        ->check(CLI::ExistingFile);

    // std::string revision = "";
    benchmark_sub->add_option("-r,--revision", revision, "Revision of this software (unique id, e.g. git commit hash)")
        ->default_val("undefined");

    // std::string machine_id = "";
    benchmark_sub->add_option("-m,--machine", machine_id, "Identifier of this computer (e.g. hostname)")
        ->default_val("undefined");

    size_t num_iterations = 1;
    benchmark_sub->add_option("-n,--iterations", num_iterations, "The number of times to run each benchmark")
        ->check(CLI::PositiveNumber)
        ->default_val(1);

    uint64_t num_warmup_iterations = 1;
    benchmark_sub
        ->add_option(
            "-w,--warmup-iterations",
            num_warmup_iterations,
            "The number of warmup iterations to run before the actual benchmarking starts"
        )
        ->check(CLI::NonNegativeNumber)
        ->default_val(1);

    benchmark_sub->callback(
        [&trees_file, &forest_file, &num_iterations, &num_warmup_iterations, &setup_results_printer]() {
            auto results_printer = setup_results_printer();

            for (uint8_t _iteration = 0; _iteration < num_iterations + num_warmup_iterations; ++_iteration) {
                bool const    warmup    = _iteration < num_warmup_iterations;
                uint8_t const iteration = warmup ? _iteration : _iteration - num_warmup_iterations;

                if (warmup) {
                    std::cerr << "[Benchmark] Running warmup iteration " << static_cast<int>(iteration) << " for file "
                              << trees_file << std::endl;
                } else {
                    std::cerr << "[Benchmark] Running iteration " << static_cast<int>(iteration) << " on file "
                              << trees_file << std::endl;
                }
                benchmark(warmup, iteration, trees_file, forest_file, results_printer);
            }
        }
    );

    // Require exactly one subcommand
    app.require_subcommand(1);

    // Parse options and run subcommand
    CLI11_PARSE(app, argc, argv);
}
