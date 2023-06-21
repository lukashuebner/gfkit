#include <fstream>
#include <sstream>
#include <string>

#include <tdt/load/forest-compressor.hpp>
#include <tdt/sequence-forest.hpp>
#include <tdt/tskit.hpp>

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"
#include "catch2/catch_approx.hpp"
#include "perf.hpp"
#include "tdt/graph/compressed-forest.hpp"
#include "tdt/load/compressed-forest-serialization.hpp"
#include "tdt/sequence/allele-frequency-spectrum.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"
#include "timer.hpp"

void benchmark(
    bool const                 warmup,
    uint8_t const              iteration,
    std::string const&         ts_file,
    ResultsPrinter&            results_printer,
    std::optional<std::string> sf_input_file,
    std::optional<std::string> sf_output_file
) {
    auto log_time = [&results_printer, iteration, &ts_file](
                        bool const         warmup,
                        std::string const& section,
                        std::string const& variant,
                        Timer::duration    duration
                    ) {
        results_printer.print_timer(warmup, section, variant, ts_file, duration, iteration);
    };

    auto log_mem = [&results_printer, iteration, &ts_file](
                       bool const                 warmup,
                       std::string const&         section,
                       std::string const&         variant,
                       MemoryUsage::Report const& report
                   ) {
        results_printer.print_memory(warmup, section, variant, ts_file, report, iteration);
    };

    // Benchmark tree sequence loading
    Timer             timer;
    MemoryUsage       memory_usage;
    TSKitTreeSequence tree_sequence(ts_file);
    do_not_optimize(tree_sequence);
    if (!warmup) {
        log_time(warmup, "load", "tskit", timer.stop());
    }

    // Benchmark building the DAG from the tree sequence
    CompressedForest       compressed_forest;
    GenomicSequenceStorage sequence_store;
    if (sf_input_file) {
        memory_usage.start();
        timer.start();

        CompressedForestIO::load(*sf_input_file, compressed_forest, sequence_store);

        log_time(warmup, "load_forest_file", "sfkit", timer.stop());
        log_mem(warmup, "load_forest_file", "sfkit", memory_usage.stop());
    } else {
        memory_usage.start();
        timer.start();

        ForestCompressor forest_compressor(tree_sequence);
        compressed_forest = forest_compressor.compress();

        log_time(warmup, "compress_forest", "sfkit", timer.stop());
        log_mem(warmup, "compress_forest", "sfkit", memory_usage.stop());

        memory_usage.start();
        timer.start();

        sequence_store = GenomicSequenceStorage(tree_sequence, forest_compressor);

        log_time(warmup, "init_genomic_sequence_storage", "sfkit", timer.stop());
        log_mem(warmup, "init_genomic_sequence_storage", "sfkit", memory_usage.stop());
    }

    if (sf_output_file) {
        memory_usage.start();
        timer.start();

        CompressedForestIO::save(*sf_output_file, compressed_forest, sequence_store);

        log_time(warmup, "save_forest_file", "sfkit", timer.stop());
        log_mem(warmup, "save_forest_file_delta_rss", "sfkit", memory_usage.stop());
    }

    // Benchmark computing subtree sizes
    memory_usage.start();
    timer.start();

    EdgeListGraph const& dag = compressed_forest.postorder_edges();
    compressed_forest.compute_num_samples_below();
    do_not_optimize(compressed_forest.num_samples_below());

    log_time(warmup, "compute_subtree_sizes", "sfkit", timer.stop());
    log_mem(warmup, "compute_subtree_sizes", "sfkit", memory_usage.stop());

    // Benchmark computing the AFS
    memory_usage.start();
    timer.start();

    AlleleFrequencySpectrum<PerfectDNAHasher> afs(sequence_store, compressed_forest);
    do_not_optimize(afs);

    log_time(warmup, "compute_afs", "sfkit", timer.stop());
    log_mem(warmup, "compute_afs", "sfkit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    auto reference_afs = tree_sequence.allele_frequency_spectrum();
    do_not_optimize(reference_afs);

    log_time(warmup, "compute_afs", "tskit", timer.stop());
    log_mem(warmup, "compute_afs", "tskit", memory_usage.stop());

    // Does our AFS match the one computed by tskit?
    bool const equal = std::ranges::equal(
        std::span(afs).subspan(1, afs.num_samples() - 1),
        std::span(reference_afs).subspan(1, afs.num_samples() - 1)
    );
    if (!equal) {
        std::cerr << "ERROR !! AFS mismatch between tskit and sf" << std::endl;
        std::exit(1);
    }

    // Benchmark computing the diversity
    SequenceForest sequence_forest(compressed_forest, sequence_store);
    memory_usage.start();
    timer.start();

    double const sfkit_diversity = sequence_forest.diversity();
    do_not_optimize(sfkit_diversity);

    log_time(warmup, "diversity", "sfkit", timer.stop());
    log_mem(warmup, "diversity", "sfkit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_diversity = tree_sequence.diversity();
    do_not_optimize(tskit_diversity);

    log_time(warmup, "diversity", "tskit", timer.stop());
    log_mem(warmup, "diversity", "tskit", memory_usage.stop());

    // Does our AFS match the one computed by tskit?
    if (sfkit_diversity == Catch::Approx(tskit_diversity).epsilon(1e-6)) {
        std::cerr << "ERROR !! Diversity mismatch between tskit and sfkit" << std::endl;
        std::cerr << "    " << tskit_diversity << " vs. " << sfkit_diversity << " (tskit vs. sfkit)" << std::endl;
        std::exit(1);
    }
}

int main(int argc, char** argv) {
    CLI::App app{"Forest Compression Benchmark"};

    size_t num_iterations = 1;
    app.add_option("-n,--iterations", num_iterations, "The number of times to run each benchmark")
        ->check(CLI::PositiveNumber)
        ->default_val(1);

    std::string ts_file = "";
    app.add_option("file,-f,--file", ts_file, "The tree sequence file to benchmark on. If not specified.")
        ->check(CLI::ExistingFile)
        ->required();

    std::string revision = "";
    app.add_option("-r,--revision", revision, "Revision of the benchmarked program (unique id, e.g. git commit hash)")
        ->default_val("undefined");

    std::string machine_id = "";
    app.add_option("-m,--machine", machine_id, "Identifier of this computer (e.g. hostname)")->default_val("undefined");

    std::optional<std::string> sf_output_file = std::nullopt;
    auto sf_output_file_opt = app.add_option("-o,--output", sf_output_file, "Output file for the compressed forest");

    std::optional<std::string> sf_input_file = std::nullopt;
    auto sf_input_file_opt = app.add_option("-i,--input", sf_input_file, "Input file for the compressed forest")
                                 ->check(CLI::ExistingFile)
                                 ->excludes(sf_output_file_opt);

    uint64_t num_warmup_iterations = 1;
    app.add_option(
           "-w,--warmup-iterations",
           num_warmup_iterations,
           "The number of warmup iterations to run before the actual benchmarking starts"
    )
        ->check(CLI::NonNegativeNumber)
        ->default_val(1);

    CLI11_PARSE(app, argc, argv);
    KASSERT(!(sf_output_file && sf_input_file), "Cannot specify both input and output file", tdt::assert::light);

    // Set-up results printer
    ResultsPrinter results_printer(std::cout, revision, machine_id);
    results_printer.print_header();

    for (uint8_t _iteration = 0; _iteration < num_iterations + num_warmup_iterations; ++_iteration) {
        bool const    warmup    = _iteration < num_warmup_iterations;
        uint8_t const iteration = warmup ? _iteration : _iteration - num_warmup_iterations;

        if (warmup) {
            std::cerr << "Running warmup iteration " << static_cast<int>(iteration) << " for file " << ts_file << "\n";
        } else {
            std::cerr << "Running iteration " << static_cast<int>(iteration) << " on file " << ts_file << "\n";
        }
        benchmark(warmup, iteration, ts_file, results_printer, sf_input_file, sf_output_file);
    }
}
