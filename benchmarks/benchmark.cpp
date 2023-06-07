#include <tdt/load/forest-compressor.hpp>
#include <tdt/tskit.hpp>

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"
#include "perf.hpp"
#include "tdt/sequence/allele-frequency-spectrum.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"

// TODO Switch to Celero or my own benchmarking framework
class Timer {
public:
    Timer() {
        start();
    }

    uint64_t stop() {
        asm("" ::: "memory");
        auto end = std::chrono::steady_clock::now();
        asm("" ::: "memory");

        return duration_cast<std::chrono::nanoseconds>(end - _start).count();
    }

    void start() {
        asm("" ::: "memory"); // prevent compiler reordering
        _start = std::chrono::steady_clock::now();
        asm("" ::: "memory");
    }

    template <class Func>
    static uint64_t time_func(Func func) {
        Timer timer;
        func();
        return timer.stop();
    }

private:
    std::chrono::time_point<std::chrono::steady_clock> _start;
};

// This assembly magic is from Google Benchmark
// see https://github.com/google/benchmark/blob/v1.7.1/include/benchmark/benchmark.h#L443-L446
template <typename T>
void do_not_optimize(T const& val) {
    asm volatile("" : : "r,m"(val) : "memory");
}

template <typename T>
void do_not_optimize(T& val) {
#if defined(__clang__)
    asm volatile("" : "+r,m"(val) : : "memory");
#else
    asm volatile("" : "+m,r"(val) : : "memory");
#endif
}

class ResultsPrinter {
public:
    ResultsPrinter(std::ostream& stream, std::string const& revision, std::string const& machine_id)
        : _stream(stream),
          _revision(revision),
          _machine_id(machine_id) {}

    void print_header() {
        _stream << "algorithm,variant,dataset,revision,machine_id,iteration,walltime_ns\n";
    }

    void print(
        bool const         warmup,
        std::string const& algorithm,
        std::string const& variant,
        std::string const& dataset,
        uint64_t           time_ns,
        size_t             iteration
    ) {
        if (!warmup) {
            _stream << algorithm << ",";
            _stream << variant << ",";
            _stream << dataset << ",";
            _stream << _revision << ",";
            _stream << _machine_id << ",";
            _stream << iteration << ",";
            _stream << time_ns << "\n";
        }
    }

private:
    std::ostream& _stream;
    std::string   _revision;
    std::string   _machine_id;
};

void benchmark(
    bool const warmup, uint8_t const iteration, std::string const& ts_file, ResultsPrinter& results_printer
) {
    // Benchmark tree sequence loading
    Timer             timer;
    TSKitTreeSequence tree_sequence(ts_file);
    do_not_optimize(tree_sequence);
    if (!warmup) {
        results_printer.print(warmup, "load", "tskit", ts_file, timer.stop(), iteration);
    }

    // Benchmark building the DAG from the tree sequence
    timer.start();
    CompressedForest       compressed_forest(tree_sequence);
    GenomicSequenceStorage sequence_store(tree_sequence, compressed_forest);
    if (!warmup) {
        results_printer.print(warmup, "build_sf", "sf", ts_file, timer.stop(), iteration);
    }

    // Benchmark computing subtree sizes
    // Compute the number of samples below each root.
    EdgeListGraph const& dag = compressed_forest.postorder_edges();
    compressed_forest.compute_num_samples_below();
    do_not_optimize(compressed_forest.num_samples_below());
    results_printer.print(warmup, "compute_subtree_sizes", "sf", ts_file, timer.stop(), iteration);

    // Benchmark computing the AFS
    // Compute the number of samples below each root.
    timer.start();
    AlleleFrequencySpectrum<PerfectDNAHasher> afs(sequence_store, compressed_forest);
    do_not_optimize(afs);
    results_printer.print(warmup, "compute_afs", "sf", ts_file, timer.stop(), iteration);

    timer.start();
    auto reference_afs = tree_sequence.allele_frequency_spectrum();
    do_not_optimize(reference_afs);
    results_printer.print(warmup, "compute_afs", "tskit", ts_file, timer.stop(), iteration);
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

    uint64_t num_warmup_iterations = 1;
    app.add_option(
           "-w,--warmup-iterations",
           num_warmup_iterations,
           "The number of warmup iterations to run before the actual benchmarking starts"
       )
        ->check(CLI::NonNegativeNumber)
        ->default_val(1);

    CLI11_PARSE(app, argc, argv);

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
        benchmark(warmup, iteration, ts_file, results_printer);
    }
}
