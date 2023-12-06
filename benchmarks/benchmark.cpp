#include <fstream>
#include <sstream>
#include <string>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <catch2/catch_approx.hpp>

#include "BenchmarkRunner.hpp"
#include "ResultsPrinter.hpp"
#include "perf.hpp"
#include "print_ds_stats.hpp"
#include "sfkit/SuccinctForest.hpp"
#include "sfkit/bp/BPForestCompressor.hpp"
#include "sfkit/dag/DAGCompressedForest.hpp"
#include "sfkit/dag/DAGForestCompressor.hpp"
#include "sfkit/io/CompressedForestIO.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"
#include "sfkit/stats/AlleleFrequencySpectrum.hpp"
#include "sfkit/tskit/tskit.hpp"
#include "timer.hpp"

void benchmark(
    bool const         warmup,
    uint8_t const      iteration,
    std::string const& trees_file,
    std::string const& dag_forest_file,
    std::string const& bp_forest_file,
    ResultsPrinter&    results_printer
) {
    auto log_time = [&results_printer,
                     iteration,
                     &trees_file,
                     warmup](std::string const& section, std::string const& variant, Timer::duration duration) {
        results_printer.print_timer(warmup, section, variant, trees_file, duration, iteration);
    };

    auto log_mem = [&results_printer,
                    iteration,
                    &trees_file,
                    warmup](std::string const& section, std::string const& variant, MemoryUsage::Report const& report) {
        results_printer.print_memory(warmup, section, variant, trees_file, report, iteration);
    };

    Timer           timer;
    MemoryUsage     memory_usage;
    BenchmarkRunner bench(memory_usage, timer, log_mem, log_time, []() { std::exit(EXIT_FAILURE); });

    // --- Benchmark tree sequence loading ---
    bench.start();
    sfkit::tskit::TSKitTreeSequence tree_sequence(trees_file);
    do_not_optimize(tree_sequence);
    bench.stop("load_trees_file", "tskit");

    // --- Benchmark loading the DAG-based tree sequence ---
    bench.start();
    sfkit::dag::DAGCompressedForest         dag_forest;
    sfkit::sequence::GenomicSequenceFactory dag_sequence_factory(tree_sequence);
    sfkit::sequence::GenomicSequence        dag_sequence;
    // TODO Write IO functionality for the SequenceForest class
    sfkit::io::CompressedForestIO::load(dag_forest_file, dag_forest, dag_sequence);
    bench.stop("load_forest_file", "sfkit_dag");

    // --- Benchmark loading the DAG-based tree sequence ---
    bench.start();
    sfkit::bp::BPCompressedForest           bp_forest;
    sfkit::sequence::GenomicSequenceFactory bp_sequence_factory(tree_sequence);
    sfkit::sequence::GenomicSequence        bp_sequence;
    sfkit::io::CompressedForestIO::load(bp_forest_file, bp_forest, bp_sequence);
    bench.stop("load_forest_file", "sfkit_bp");

    // --- Benchmark computing subtree sizes ---
    // TODO Add custom checker? Overload operator==? (Approx checking is a little bit more involved)
    // TODO Add the numsamples below function to SuccinctForest

    using SetOfSampleSets = sfkit::sequence::SetOfSampleSets<1>;
    auto all_samples      = dag_forest.all_samples();

    bench.start();
    auto num_samples_below = sfkit::sequence::NumSamplesBelowFactory::build(dag_forest, all_samples);
    do_not_optimize(num_samples_below);
    bench.stop("compute_subtree_sizes", "sfkit_dag");

    bench.start();
    auto bp_num_samples_below = sfkit::sequence::NumSamplesBelowFactory::build(bp_forest, all_samples);
    do_not_optimize(bp_num_samples_below);
    bench.stop("compute_subtree_sizes", "sfkit_bp");

    // --- Construct the SequenceForest object on which we then call the high-level operations. ---
    sfkit::DAGSuccinctForest dag_succinct_forest(std::move(dag_forest), std::move(dag_sequence));

    sfkit::BPSuccinctForest bp_succinct_forest(std::move(bp_forest), std::move(bp_sequence));

    // --- Output some statistics ---
    print_ds_stats(
        warmup,
        trees_file,
        iteration,
        results_printer,
        dag_succinct_forest,
        bp_succinct_forest,
        tree_sequence
    );

    // --- Benchmark computing the AFS ---
    bench.start();
    auto const sfkit_afs = dag_succinct_forest.allele_frequency_spectrum();
    do_not_optimize(sfkit_afs);
    bench.stop("afs", "sfkit_dag");

    bench.start();
    auto const tskit_afs = tree_sequence.allele_frequency_spectrum();
    do_not_optimize(tskit_afs);
    bench.stop("afs", "tskit");

    bench.start();
    auto const sfkit_bp_afs = bp_succinct_forest.allele_frequency_spectrum();
    do_not_optimize(sfkit_bp_afs);
    bench.stop("afs", "sfkit_bp");

    // Does our AFS match the one computed by tskit?
    for (size_t i = 1; i < sfkit_afs.num_samples() - 1; ++i) {
        if (sfkit_afs[i] != tskit_afs[i]) {
            std::cerr << "ERROR !! AFS mismatch between tskit and sfkit" << std::endl;
            std::cerr << "    " << tskit_afs[i] << " vs. " << sfkit_afs[i] << " (tskit vs. sfkit)" << std::endl;
            std::exit(1);
        }

        if (sfkit_bp_afs[i] != tskit_afs[i]) {
            std::cerr << "ERROR !! AFS mismatch between tskit and sfkit_bp" << std::endl;
            std::cerr << "    " << tskit_afs[i] << " vs. " << sfkit_bp_afs[i] << " (tskit vs. sfkit_bp)" << std::endl;
            std::exit(1);
        }
    }

    // --- Benchmark computing the divergence ---
    // TODO Extract helpers for sample set creation (also used in unit tests)
    sfkit::samples::SampleSet sample_set_1(dag_succinct_forest.num_samples());
    sfkit::samples::SampleSet sample_set_2(dag_succinct_forest.num_samples());
    bool                      flip = false;
    for (sfkit::samples::SampleId sample: dag_succinct_forest.all_samples()) {
        if (flip) {
            sample_set_1.add(sample);
        } else {
            sample_set_2.add(sample);
        }
        flip = !flip;
    }

    bench.bench_and_check(
        "divergence",
        [&tree_sequence, &sample_set_1, &sample_set_2]() {
            return tree_sequence.divergence(sample_set_1, sample_set_2);
        },
        [&dag_succinct_forest, &sample_set_1, &sample_set_2]() {
            return dag_succinct_forest.divergence(sample_set_1, sample_set_2);
        },
        [&bp_succinct_forest, &sample_set_1, &sample_set_2]() {
            return bp_succinct_forest.divergence(sample_set_1, sample_set_2);
        }
    );

    // --- Benchmark computing Patterson's f{2,3,4} ---
    constexpr size_t                       num_sample_sets = 4;
    std::vector<sfkit::samples::SampleSet> sample_sets(num_sample_sets, dag_succinct_forest.num_samples());

    size_t idx = 0;
    for (sfkit::samples::SampleId sample: dag_succinct_forest.forest().leaves()) {
        KASSERT(sample < dag_succinct_forest.num_samples(), "Sample id out of range", sfkit::assert::light);
        sample_sets[idx].add(sample);
        idx = (idx + 1ul) % num_sample_sets;
    }

    bench.bench_and_check(
        "f4",
        [&tree_sequence, &sample_sets]() {
            return tree_sequence.f4(sample_sets[0], sample_sets[1], sample_sets[2], sample_sets[3]);
        },
        [&dag_succinct_forest, &sample_sets]() {
            return dag_succinct_forest.f4(sample_sets[0], sample_sets[1], sample_sets[2], sample_sets[3]);
        },
        [&bp_succinct_forest, &sample_sets]() {
            return bp_succinct_forest.f4(sample_sets[0], sample_sets[1], sample_sets[2], sample_sets[3]);
        }
    );

    bench.bench_and_check(
        "f3",
        [&tree_sequence, &sample_sets]() { return tree_sequence.f3(sample_sets[0], sample_sets[1], sample_sets[2]); },
        [&dag_succinct_forest, &sample_sets]() {
            return dag_succinct_forest.f3(sample_sets[0], sample_sets[1], sample_sets[2]);
        },
        [&bp_succinct_forest, &sample_sets]() {
            return bp_succinct_forest.f3(sample_sets[0], sample_sets[1], sample_sets[2]);
        }
    );

    bench.bench_and_check(
        "f2",
        [&tree_sequence, &sample_sets]() { return tree_sequence.f2(sample_sets[0], sample_sets[1]); },
        [&dag_succinct_forest, &sample_sets]() { return dag_succinct_forest.f2(sample_sets[0], sample_sets[1]); },
        [&bp_succinct_forest, &sample_sets]() { return bp_succinct_forest.f2(sample_sets[0], sample_sets[1]); }
    );

    // --- Benchmark computing the diversity ---
    bench.bench_and_check(
        "diversity",
        [&tree_sequence]() { return tree_sequence.diversity(); },
        [&dag_succinct_forest]() { return dag_succinct_forest.diversity(); },
        [&bp_succinct_forest]() { return bp_succinct_forest.diversity(); }
    );

    // --- Benchmark computing the number of segregating sites ---
    bench.bench_and_check(
        "num_segregating_sites",
        [&tree_sequence]() { return tree_sequence.num_segregating_sites(); },
        [&dag_succinct_forest]() { return dag_succinct_forest.num_segregating_sites(); },
        [&bp_succinct_forest]() { return bp_succinct_forest.num_segregating_sites(); }
    );

    // --- Benchmark computing Tajima's D. ---
    // tskit does implement this only in Python, not in C. We're using
    // experiments/scripts/benchmark-tskits-tajimas-d.py to measure tskit's performance.
    bench.start();
    double const tajimas_d = dag_succinct_forest.tajimas_d();
    do_not_optimize(tajimas_d);
    bench.stop("tajimas_d", "sfkit_dag");

    results_printer.print(warmup, "tajimas_d", "sfkit_dag", trees_file, "tajimas_d", tajimas_d, "1", iteration);

    // --- Benchmark computing the FST. ---
    // tskit's implementation is in Python; see Tajima's D above.
    bench.start();
    double const sfkit_fst = dag_succinct_forest.fst(sample_set_1, sample_set_2);
    do_not_optimize(sfkit_fst);
    bench.stop("fst", "sfkit_dag");

    results_printer.print(warmup, "fst", "sfkit_dag", trees_file, "fst", sfkit_fst, "1", iteration);
}
