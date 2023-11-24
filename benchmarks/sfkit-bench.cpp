#include <fstream>
#include <sstream>
#include <string>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <catch2/catch_approx.hpp>

#include "perf.hpp"
#include "sfkit/SuccinctForest.hpp"
#include "sfkit/graph/CompressedForest.hpp"
#include "sfkit/load/BPForestCompressor.hpp"
#include "sfkit/load/CompressedForestIO.hpp"
#include "sfkit/load/ForestCompressor.hpp"
#include "sfkit/sequence/AlleleFrequencySpectrum.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"
#include "sfkit/tskit.hpp"
#include "timer.hpp"

constexpr double FLOAT_EQ_EPS = 1e-4;

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

    sfkit::io::CompressedForestIO::save(forest_file, forest, sequence);

    log_time("save_forest_file", "sfkit", timer.stop());
    log_mem("save_forest_file", "sfkit", memory_usage.stop());

    // Compress the .trees to .bpforest file

    memory_usage.start();
    timer.start();

    using SetOfSampleSets = NumSamplesBelow<1>::SetOfSampleSets;
    BPForestCompressor     bp_forest_compressor(tree_sequence);
    GenomicSequenceFactory bp_sequence_factory(tree_sequence);
    BPCompressedForest     bp_forest   = bp_forest_compressor.compress(bp_sequence_factory);
    GenomicSequence        bp_sequence = bp_sequence_factory.move_storage();

    log_time("compress_forest_and_sequence", "bp-sfkit", timer.stop());
    log_mem("compress_forest_and_sequence", "bp-sfkit", memory_usage.stop());

    // Save the BP-compressed forest and sequence to a .bpforest file
    memory_usage.start();
    timer.start();

    sfkit::io::CompressedForestIO::save(bp_forest_file, bp_forest, bp_sequence);

    log_time("save_forest_file", "bp-sfkit", timer.stop());
    log_mem("save_forest_file", "bp-sfkit", memory_usage.stop());
}

void benchmark(
    bool const         warmup,
    uint8_t const      iteration,
    std::string const& trees_file,
    std::string const& forest_file,
    std::string const& bp_forest_file,
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

    // --- Benchmark tree sequence loading ---
    Timer             timer;
    MemoryUsage       memory_usage;
    TSKitTreeSequence tree_sequence(trees_file);
    do_not_optimize(tree_sequence);
    if (!warmup) {
        log_time(warmup, "load_trees_file", "tskit", timer.stop());
    }

    // --- Benchmark loading the DAG-based tree sequence ---
    memory_usage.start();
    timer.start();

    CompressedForest       forest;
    GenomicSequenceFactory sequence_factory(tree_sequence);
    GenomicSequence        sequence;

    // TODO Write IO functionality for the SequenceForest class
    sfkit::io::CompressedForestIO::load(forest_file, forest, sequence);

    log_time(warmup, "load_forest_file", "sfkit", timer.stop());
    log_mem(warmup, "load_forest_file", "sfkit", memory_usage.stop());

    // --- Benchmark loading the DAG-based tree sequence ---
    memory_usage.start();
    timer.start();

    BPCompressedForest     bp_forest;
    GenomicSequenceFactory bp_sequence_factory(tree_sequence);
    GenomicSequence        bp_sequence;

    sfkit::io::CompressedForestIO::load(bp_forest_file, bp_forest, bp_sequence);

    log_time(warmup, "load_forest_file", "bp-sfkit", timer.stop());
    log_mem(warmup, "load_forest_file", "bp-sfkit", memory_usage.stop());

    // --- Benchmark computing subtree sizes ---
    // memory_usage.start();
    // timer.start();

    // EdgeListGraph const& dag = forest.postorder_edges();
    // using SetOfSampleSets    = NumSamplesBelow<1>::SetOfSampleSets;
    auto all_samples = forest.all_samples();
    // auto num_samples_below   = NumSamplesBelow<1>(dag, SetOfSampleSets{std::cref(all_samples)});
    // do_not_optimize(num_samples_below);

    // log_time(warmup, "compute_subtree_sizes", "sfkit", timer.stop());
    // log_mem(warmup, "compute_subtree_sizes", "sfkit", memory_usage.stop());

    // --- Benchmark computing subtree sizes using the BP-based data structure---
    memory_usage.start();
    timer.start();

    for (int i = 0; i < 100; i++) {
        auto bp_num_samples_below = NumSamplesBelowFactory::build(bp_forest, all_samples);
        do_not_optimize(bp_num_samples_below);
    }
    return;

    log_time(warmup, "compute_subtree_sizes", "sfkit-bp", timer.stop());
    log_mem(warmup, "compute_subtree_sizes", "sfkit-bp", memory_usage.stop());

    // --- Construct the SequenceForest object on which we then call the high-level operations. ---
    // TODO The SequenceForest does not need to hold the tree_sequence at all times; it's only used during
    // construction. (Check this!)
    SuccinctForest sequence_forest(std::move(tree_sequence), std::move(forest), std::move(sequence));

    // --- Output some statistics ---
    results_printer
        .print(warmup, "stats", "sfkit", trees_file, "num_trees", sequence_forest.num_trees(), "1", iteration);
    results_printer.print(
        warmup,
        "stats",
        "tskit",
        trees_file,
        "num_trees",
        sequence_forest.tree_sequence().num_trees(),
        "1",
        iteration
    );

    results_printer
        .print(warmup, "stats", "sfkit", trees_file, "num_mutations", sequence_forest.num_mutations(), "1", iteration);
    results_printer.print(
        warmup,
        "stats",
        "tskit",
        trees_file,
        "num_mutations",
        sequence_forest.tree_sequence().num_mutations(),
        "1",
        iteration
    );

    results_printer
        .print(warmup, "stats", "sfkit", trees_file, "num_sites", sequence_forest.num_sites(), "1", iteration);
    results_printer.print(
        warmup,
        "stats",
        "tskit",
        trees_file,
        "num_sites",
        sequence_forest.tree_sequence().num_sites(),
        "1",
        iteration
    );

    results_printer
        .print(warmup, "stats", "sfkit", trees_file, "num_samples", sequence_forest.num_samples(), "1", iteration);
    results_printer.print(
        warmup,
        "stats",
        "tskit",
        trees_file,
        "num_samples",
        sequence_forest.tree_sequence().num_samples(),
        "1",
        iteration
    );

    results_printer.print(
        warmup,
        "stats",
        "sfkit",
        trees_file,
        "num_unique_subtrees",
        sequence_forest.num_unique_subtrees(),
        "1",
        iteration
    );

    results_printer.print(
        warmup,
        "stats",
        "sfkit",
        trees_file,
        "num_subtrees_with_mutations",
        sequence_forest.num_subtrees_with_mutations(),
        "1",
        iteration
    );

    results_printer.print(
        warmup,
        "stats",
        "tskit",
        trees_file,
        "num_edges",
        sequence_forest.tree_sequence().num_edges(),
        "1",
        iteration
    );

    // --- Benchmark computing the AFS ---
    memory_usage.start();
    timer.start();

    auto const sfkit_afs = sequence_forest.allele_frequency_spectrum();
    do_not_optimize(sfkit_afs);

    log_time(warmup, "afs", "sfkit", timer.stop());
    log_mem(warmup, "afs", "sfkit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    auto const tskit_afs = sequence_forest.tree_sequence().allele_frequency_spectrum();
    do_not_optimize(tskit_afs);

    log_time(warmup, "afs", "tskit", timer.stop());
    log_mem(warmup, "afs", "tskit", memory_usage.stop());

    // Does our AFS match the one computed by tskit?
    for (size_t i = 1; i < sfkit_afs.num_samples() - 1; ++i) {
        if (sfkit_afs[i] != tskit_afs[i]) {
            std::cerr << "ERROR !! AFS mismatch between tskit and sfkit" << std::endl;
            std::cerr << "    " << tskit_afs[i] << " vs. " << sfkit_afs[i] << " (tskit vs. sfkit)" << std::endl;
            std::exit(1);
        }
    }

    // --- Benchmark computing the divergence ---
    SampleSet sample_set_1(sequence_forest.num_samples());
    SampleSet sample_set_2(sequence_forest.num_samples());
    bool      flip = false;
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

    // --- Benchmark computing Patterson's f{2,3,4} ---
    constexpr size_t       num_sample_sets = 4;
    std::vector<SampleSet> sample_sets(num_sample_sets, sequence_forest.num_samples());

    size_t idx = 0;
    for (SampleId sample: sequence_forest.forest().leaves()) {
        KASSERT(sample < sequence_forest.num_samples(), "Sample id out of range", sfkit::assert::light);
        sample_sets[idx].add(sample);
        idx = (idx + 1ul) % num_sample_sets;
    }

    memory_usage.start();
    timer.start();

    double const sfkit_f4 = sequence_forest.f4(sample_sets[0], sample_sets[1], sample_sets[2], sample_sets[3]);

    log_time(warmup, "f4", "sfkit", timer.stop());
    log_mem(warmup, "f4", "sfkit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_f4 =
        sequence_forest.tree_sequence().f4(sample_sets[0], sample_sets[1], sample_sets[2], sample_sets[3]);

    log_time(warmup, "f4", "tskit", timer.stop());
    log_mem(warmup, "f4", "tskit", memory_usage.stop());

    if (sfkit_f4 != Catch::Approx(tskit_f4).epsilon(1e-4)) {
        std::cerr << "ERROR !! f4 mismatch between tskit and sfkit" << std::endl;
        std::cerr << "    " << tskit_f4 << " vs. " << sfkit_f4 << " (tskit vs. sfkit)" << std::endl;
        // std::exit(1);
    }

    memory_usage.start();
    timer.start();

    double const sfkit_f3 = sequence_forest.f3(sample_sets[0], sample_sets[1], sample_sets[2]);

    log_time(warmup, "f3", "sfkit", timer.stop());
    log_mem(warmup, "f3", "sfkit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_f3 = sequence_forest.tree_sequence().f3(sample_sets[0], sample_sets[1], sample_sets[2]);

    log_time(warmup, "f3", "tskit", timer.stop());
    log_mem(warmup, "f3", "tskit", memory_usage.stop());

    if (sfkit_f3 != Catch::Approx(tskit_f3).epsilon(1e-4)) {
        std::cerr << "ERROR !! f3 mismatch between tskit and sfkit" << std::endl;
        std::cerr << "    " << tskit_f3 << " vs. " << sfkit_f3 << " (tskit vs. sfkit)" << std::endl;
        // std::exit(1);
    }

    memory_usage.start();
    timer.start();

    double const sfkit_f2 = sequence_forest.f2(sample_sets[0], sample_sets[1]);

    log_time(warmup, "f2", "sfkit", timer.stop());
    log_mem(warmup, "f2", "sfkit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_f2 = sequence_forest.tree_sequence().f2(sample_sets[0], sample_sets[1]);

    log_time(warmup, "f2", "tskit", timer.stop());
    log_mem(warmup, "f2", "tskit", memory_usage.stop());

    if (sfkit_f2 != Catch::Approx(tskit_f2).epsilon(1e-4)) {
        std::cerr << "ERROR !! f2 mismatch between tskit and sfkit" << std::endl;
        std::cerr << "    " << tskit_f2 << " vs. " << sfkit_f2 << " (tskit vs. sfkit)" << std::endl;
        // std::exit(1);
    }

    // --- Benchmark computing the diversity ---
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

    // --- Benchmark computing the number of segregating sites ---
    memory_usage.start();
    timer.start();

    uint64_t const sfkit_num_seg_sites = sequence_forest.num_segregating_sites();
    do_not_optimize(sfkit_num_seg_sites);

    log_time(warmup, "num_segregating_sites", "sfkit", timer.stop());
    log_mem(warmup, "num_segregating_sites", "sfkit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_num_seg_sites = sequence_forest.tree_sequence().num_segregating_sites();
    do_not_optimize(tskit_num_seg_sites);

    log_time(warmup, "num_segregating_sites", "tskit", timer.stop());
    log_mem(warmup, "num_segregating_sites", "tskit", memory_usage.stop());

    // Does our number of segregating sites match the one computed by tskit?
    if (sfkit_num_seg_sites != Catch::Approx(tskit_num_seg_sites).epsilon(FLOAT_EQ_EPS)) {
        std::cerr << "ERROR !! Number of segregating sites mismatch between tskit and sfkit" << std::endl;
        std::cerr << "    " << tskit_num_seg_sites << " vs. " << sfkit_num_seg_sites << " (tskit vs. sfkit)"
                  << std::endl;
        std::exit(1);
    }

    // --- Benchmark computing Tajima's D. ---
    // tskit does implement this only in Python, not in C. We're using
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

    // --- Benchmark computing the FST. ---
    // tskit's implementation is in Python; see Tajima's D above.
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

    std::string bp_forest_file = "";
    compress_sub->add_option("-b,--bp-forest-file", bp_forest_file, "Input file for the BP-compressed forest")
        ->check(CLI::NonexistentPath);

    compress_sub->add_option("-r,--revision", revision, "Revision of this software (unique id, e.g. git commit hash)")
        ->default_val("undefined");

    compress_sub->add_option("-m,--machine", machine_id, "Identifier of this computer (e.g. hostname)")
        ->default_val("undefined");

    compress_sub->callback([&trees_file, &forest_file, &bp_forest_file, &setup_results_printer]() {
        auto results_printer = setup_results_printer();
        std::cerr << "Compressing tree sequence " << trees_file << " -> " << forest_file << std::endl;
        compress(trees_file, forest_file, bp_forest_file, results_printer);
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

    // std::string bp_forest_file = "";
    benchmark_sub->add_option("-b,--bp-forest-file", bp_forest_file, "Input file for the BP-compressed forest")
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
        [&trees_file, &forest_file, &bp_forest_file, &num_iterations, &num_warmup_iterations, &setup_results_printer](
        ) {
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
                benchmark(warmup, iteration, trees_file, forest_file, bp_forest_file, results_printer);
            }
        }
    );

    // Require exactly one subcommand
    app.require_subcommand(1);

    // Parse options and run subcommand
    CLI11_PARSE(app, argc, argv);
}