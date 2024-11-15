#include <fstream>
#include <sstream>
#include <string>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <catch2/catch_approx.hpp>
#include <sfkit/SuccinctForest.hpp>
#include <sfkit/bp/BPForestCompressor.hpp>
#include <sfkit/dag/DAGCompressedForest.hpp>
#include <sfkit/dag/DAGForestCompressor.hpp>
#include <sfkit/io/CompressedForestIO.hpp>
#include <sfkit/sequence/GenomicSequence.hpp>
#include <sfkit/stats/AlleleFrequencySpectrum.hpp>
#include <sfkit/tskit/tskit.hpp>

#include "ResultsPrinter.hpp"
#include "benchmark.hpp"
#include "compress.hpp"
#include "dataset_stats.hpp"
#include "perf.hpp"
#include "timer.hpp"

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

    // Compress subcommand
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
        if (forest_file == "" && bp_forest_file == "") {
            std::cerr << "Please provide one or both of --forest-file or --bp-forest-file" << std::endl;
            return EXIT_FAILURE;
        }

        std::cerr << "Compressing tree sequence " << trees_file << std::endl;

        auto results_printer = setup_results_printer();
        compress(trees_file, forest_file, bp_forest_file, results_printer);
    
        return EXIT_SUCCESS;
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
            if (forest_file == "" && bp_forest_file == "") {
                std::cerr << "Please provide one or both of --forest-file or --bp-forest-file" << std::endl;
                return EXIT_FAILURE;
            }

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

            return EXIT_SUCCESS;
        }
    );

    // Stats subcommand
    auto stats_sub = app.add_subcommand("stats", "Output statistics on the datasets, tskit and sfkit encodings.");

    // std::string trees_file = "";
    stats_sub->add_option("-f,--trees-file", trees_file, "The tree sequence file")
        ->check(CLI::ExistingFile)
        ->required();

    // std::string forest_file = "";
    stats_sub->add_option("-i,--forest-file", forest_file, "Compressed forest file.")->check(CLI::ExistingFile);

    // std::string bp_forest_file = "";
    stats_sub->add_option("-b,--bp-forest-file", bp_forest_file, "BP-compressed forest file")->check(CLI::ExistingFile);

    // std::string revision = "";
    stats_sub->add_option("-r,--revision", revision, "Revision of this software (unique id, e.g. git commit hash)")
        ->default_val("undefined");

    // std::string machine_id = "";
    stats_sub->add_option("-m,--machine", machine_id, "Identifier of this computer (e.g. hostname)")
        ->default_val("undefined");

    stats_sub->callback([&trees_file, &forest_file, &bp_forest_file, &setup_results_printer]() {
        auto results_printer = setup_results_printer();

        dataset_stats(trees_file, forest_file, bp_forest_file, results_printer);
    });

    // Require exactly one subcommand
    app.require_subcommand(1);

    // Parse options and run subcommand
    CLI11_PARSE(app, argc, argv);
}
