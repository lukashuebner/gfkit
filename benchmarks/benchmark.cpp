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

constexpr double FLOAT_EQ_EPS = 1e-4;

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
    Timer                           timer;
    MemoryUsage                     memory_usage;
    sfkit::tskit::TSKitTreeSequence tree_sequence(trees_file);
    do_not_optimize(tree_sequence);
    if (!warmup) {
        log_time(warmup, "load_trees_file", "tskit", timer.stop());
    }

    // --- Benchmark loading the DAG-based tree sequence ---
    memory_usage.start();
    timer.start();

    sfkit::dag::DAGCompressedForest         dag_forest;
    sfkit::sequence::GenomicSequenceFactory dag_sequence_factory(tree_sequence);
    sfkit::sequence::GenomicSequence        dag_sequence;

    // TODO Write IO functionality for the SequenceForest class
    sfkit::io::CompressedForestIO::load(forest_file, dag_forest, dag_sequence);

    log_time(warmup, "load_forest_file", "sfkit_dag", timer.stop());
    log_mem(warmup, "load_forest_file", "sfkit_dag", memory_usage.stop());

    // --- Benchmark loading the DAG-based tree sequence ---
    memory_usage.start();
    timer.start();

    sfkit::bp::BPCompressedForest           bp_forest;
    sfkit::sequence::GenomicSequenceFactory bp_sequence_factory(tree_sequence);
    sfkit::sequence::GenomicSequence        bp_sequence;

    sfkit::io::CompressedForestIO::load(bp_forest_file, bp_forest, bp_sequence);

    log_time(warmup, "load_forest_file", "sfkit_bp", timer.stop());
    log_mem(warmup, "load_forest_file", "sfkit_bp", memory_usage.stop());

    // --- Benchmark computing subtree sizes ---
    memory_usage.start();
    timer.start();

    // TODO Add the numsamples below function to SuccinctForest
    using SetOfSampleSets  = sfkit::sequence::SetOfSampleSets<1>;
    auto all_samples       = dag_forest.all_samples();
    auto num_samples_below = sfkit::sequence::NumSamplesBelowFactory::build(dag_forest, all_samples);
    do_not_optimize(num_samples_below);

    log_time(warmup, "compute_subtree_sizes", "sfkit_dag", timer.stop());
    log_mem(warmup, "compute_subtree_sizes", "sfkit_dag", memory_usage.stop());

    // TODO Output CSV line
    std::cout << "# uniq subtrees | BP: " << bp_forest.num_unique_subtrees()
              << " DAG: " << dag_forest.num_unique_subtrees() << std::endl;
    std::cout << "# backrefs | BP: " << bp_forest.num_backrefs() << std::endl;
    std::cout << "# nodes | BP: " << bp_forest.num_nodes() << " DAG: " << dag_forest.num_nodes() << std::endl;
    std::cout << "# leaves | BP: " << bp_forest.num_leaves() << " DAG: " << dag_forest.num_leaves() << std::endl;
    std::cout << "# samples | BP: " << bp_forest.num_samples() << " DAG: " << dag_forest.num_samples() << std::endl;
    std::cout << "# trees | BP: " << bp_forest.num_trees() << " DAG: " << dag_forest.num_trees()
              << " TS: " << tree_sequence.num_trees() << std::endl;

    // --- Benchmark computing subtree sizes using the BP-based data structure---

    memory_usage.start();
    timer.start();

    auto bp_num_samples_below = sfkit::sequence::NumSamplesBelowFactory::build(bp_forest, all_samples);
    do_not_optimize(bp_num_samples_below);

    log_time(warmup, "compute_subtree_sizes", "sfkit_bp", timer.stop());
    log_mem(warmup, "compute_subtree_sizes", "sfkit_bp", memory_usage.stop());

    // --- Construct the SequenceForest object on which we then call the high-level operations. ---
    // TODO The SequenceForest does not need to hold the tree_sequence at all times; it's only used during
    // construction. (Check this!)
    // TODO Create shorthands for this?
    sfkit::SuccinctForest<sfkit::dag::DAGCompressedForest, sfkit::sequence::PerfectDNAHasher> dag_succinct_forest(
        std::move(dag_forest),
        std::move(dag_sequence)
    );

    sfkit::SuccinctForest<sfkit::bp::BPCompressedForest, sfkit::sequence::PerfectDNAHasher> bp_succinct_forest(
        std::move(bp_forest),
        std::move(bp_sequence)
    );

    // --- Output some statistics ---
    results_printer
        .print(warmup, "stats", "sfkit_dag", trees_file, "num_trees", dag_succinct_forest.num_trees(), "1", iteration);
    results_printer.print(warmup, "stats", "tskit", trees_file, "num_trees", tree_sequence.num_trees(), "1", iteration);

    results_printer.print(
        warmup,
        "stats",
        "sfkit_dag",
        trees_file,
        "num_mutations",
        dag_succinct_forest.num_mutations(),
        "1",
        iteration
    );
    results_printer
        .print(warmup, "stats", "tskit", trees_file, "num_mutations", tree_sequence.num_mutations(), "1", iteration);

    results_printer
        .print(warmup, "stats", "sfkit_dag", trees_file, "num_sites", dag_succinct_forest.num_sites(), "1", iteration);
    results_printer.print(warmup, "stats", "tskit", trees_file, "num_sites", tree_sequence.num_sites(), "1", iteration);

    results_printer.print(
        warmup,
        "stats",
        "sfkit_dag",
        trees_file,
        "num_samples",
        dag_succinct_forest.num_samples(),
        "1",
        iteration
    );
    results_printer
        .print(warmup, "stats", "tskit", trees_file, "num_samples", tree_sequence.num_samples(), "1", iteration);

    results_printer.print(
        warmup,
        "stats",
        "sfkit_dag",
        trees_file,
        "num_unique_subtrees",
        dag_succinct_forest.num_unique_subtrees(),
        "1",
        iteration
    );

    results_printer.print(
        warmup,
        "stats",
        "sfkit_dag",
        trees_file,
        "num_subtrees_with_mutations",
        dag_succinct_forest.num_subtrees_with_mutations(),
        "1",
        iteration
    );

    results_printer.print(warmup, "stats", "tskit", trees_file, "num_edges", tree_sequence.num_edges(), "1", iteration);

    // --- Benchmark computing the AFS ---
    memory_usage.start();
    timer.start();

    auto const sfkit_afs = dag_succinct_forest.allele_frequency_spectrum();
    do_not_optimize(sfkit_afs);

    log_time(warmup, "afs", "sfkit_dag", timer.stop());
    log_mem(warmup, "afs", "sfkit_dag", memory_usage.stop());

    memory_usage.start();
    timer.start();

    auto const tskit_afs = tree_sequence.allele_frequency_spectrum();
    do_not_optimize(tskit_afs);

    log_time(warmup, "afs", "tskit", timer.stop());
    log_mem(warmup, "afs", "tskit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    auto const sfkit_bp_afs = bp_succinct_forest.allele_frequency_spectrum();

    log_time(warmup, "afs", "sfkit_bp", timer.stop());
    log_mem(warmup, "afs", "sfkit_bp", memory_usage.stop());

    // Does our AFS match the one computed by tskit?
    // TODO Rename sfkit -> sfkit-dag
    // TODO Abstract away comparison
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
    memory_usage.start();
    timer.start();

    double const sfkit_dag_divergence = dag_succinct_forest.divergence(sample_set_1, sample_set_2);
    do_not_optimize(sfkit_dag_divergence);

    log_time(warmup, "divergence", "sfkit_dag", timer.stop());
    log_mem(warmup, "divergence", "sfkit_dag", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_divergence = tree_sequence.divergence(sample_set_1, sample_set_2);
    do_not_optimize(tskit_divergence);

    log_time(warmup, "divergence", "tskit", timer.stop());
    log_mem(warmup, "divergence", "tskit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const sfkit_bp_divergence = dag_succinct_forest.divergence(sample_set_1, sample_set_2);
    do_not_optimize(sfkit_bp_divergence);

    log_time(warmup, "divergence", "sfkit_bp", timer.stop());
    log_mem(warmup, "divergence", "sfkit_bp", memory_usage.stop());

    // Does our divergence match the one computed by tskit?
    if (sfkit_dag_divergence != Catch::Approx(tskit_divergence).epsilon(FLOAT_EQ_EPS)) {
        std::cerr << "ERROR !! Divergence mismatch between tskit and sfkit (DAG)" << std::endl;
        std::cerr << "    " << tskit_divergence << " vs. " << sfkit_dag_divergence << " (tskit vs. sfkit)" << std::endl;
        std::exit(1);
    }

    if (sfkit_bp_divergence != Catch::Approx(tskit_divergence).epsilon(FLOAT_EQ_EPS)) {
        std::cerr << "ERROR !! Divergence mismatch between tskit and sfkit (BP)" << std::endl;
        std::cerr << "    " << tskit_divergence << " vs. " << sfkit_bp_divergence << " (tskit vs. sfkit)" << std::endl;
        std::exit(1);
    }

    // --- Benchmark computing Patterson's f{2,3,4} ---
    constexpr size_t                       num_sample_sets = 4;
    std::vector<sfkit::samples::SampleSet> sample_sets(num_sample_sets, dag_succinct_forest.num_samples());

    size_t idx = 0;
    for (sfkit::samples::SampleId sample: dag_succinct_forest.forest().leaves()) {
        KASSERT(sample < dag_succinct_forest.num_samples(), "Sample id out of range", sfkit::assert::light);
        sample_sets[idx].add(sample);
        idx = (idx + 1ul) % num_sample_sets;
    }

    memory_usage.start();
    timer.start();

    double const sfkit_dag_f4 = dag_succinct_forest.f4(sample_sets[0], sample_sets[1], sample_sets[2], sample_sets[3]);

    log_time(warmup, "f4", "sfkit_dag", timer.stop());
    log_mem(warmup, "f4", "sfkit_dag", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_f4 = tree_sequence.f4(sample_sets[0], sample_sets[1], sample_sets[2], sample_sets[3]);

    log_time(warmup, "f4", "tskit", timer.stop());
    log_mem(warmup, "f4", "tskit", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const sfkit_bp_f4 = dag_succinct_forest.f4(sample_sets[0], sample_sets[1], sample_sets[2], sample_sets[3]);

    log_time(warmup, "f4", "sfkit_bp", timer.stop());
    log_mem(warmup, "f4", "sfkit_bp", memory_usage.stop());

    if (sfkit_dag_f4 != Catch::Approx(tskit_f4).epsilon(1e-4)) {
        std::cerr << "ERROR !! f4 mismatch between tskit and sfkit (DAG)" << std::endl;
        std::cerr << "    " << tskit_f4 << " vs. " << sfkit_dag_f4 << " (tskit vs. sfkit)" << std::endl;
        std::exit(1);
    }

    if (sfkit_bp_f4 != Catch::Approx(tskit_f4).epsilon(1e-4)) {
        std::cerr << "ERROR !! f4 mismatch between tskit and sfkit (BP)" << std::endl;
        std::cerr << "    " << tskit_f4 << " vs. " << sfkit_bp_f4 << " (tskit vs. sfkit)" << std::endl;
        std::exit(1);
    }

    memory_usage.start();
    timer.start();

    double const sfkit_f3 = dag_succinct_forest.f3(sample_sets[0], sample_sets[1], sample_sets[2]);

    log_time(warmup, "f3", "sfkit_dag", timer.stop());
    log_mem(warmup, "f3", "sfkit_dag", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_f3 = tree_sequence.f3(sample_sets[0], sample_sets[1], sample_sets[2]);

    log_time(warmup, "f3", "tskit", timer.stop());
    log_mem(warmup, "f3", "tskit", memory_usage.stop());

    if (sfkit_f3 != Catch::Approx(tskit_f3).epsilon(1e-4)) {
        std::cerr << "ERROR !! f3 mismatch between tskit and sfkit" << std::endl;
        std::cerr << "    " << tskit_f3 << " vs. " << sfkit_f3 << " (tskit vs. sfkit)" << std::endl;
        // std::exit(1);
    }

    memory_usage.start();
    timer.start();

    double const sfkit_f2 = dag_succinct_forest.f2(sample_sets[0], sample_sets[1]);

    log_time(warmup, "f2", "sfkit_dag", timer.stop());
    log_mem(warmup, "f2", "sfkit_dag", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_f2 = tree_sequence.f2(sample_sets[0], sample_sets[1]);

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

    double const sfkit_diversity = dag_succinct_forest.diversity();
    do_not_optimize(sfkit_diversity);

    log_time(warmup, "diversity", "sfkit_dag", timer.stop());
    log_mem(warmup, "diversity", "sfkit_dag", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_diversity = tree_sequence.diversity();
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

    uint64_t const sfkit_num_seg_sites = dag_succinct_forest.num_segregating_sites();
    do_not_optimize(sfkit_num_seg_sites);

    log_time(warmup, "num_segregating_sites", "sfkit_dag", timer.stop());
    log_mem(warmup, "num_segregating_sites", "sfkit_dag", memory_usage.stop());

    memory_usage.start();
    timer.start();

    double const tskit_num_seg_sites = tree_sequence.num_segregating_sites();
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

    double const tajimas_d = dag_succinct_forest.tajimas_d();
    do_not_optimize(tajimas_d);

    log_time(warmup, "tajimas_d", "sfkit_dag", timer.stop());
    log_mem(warmup, "tajimas_d", "sfkit_dag", memory_usage.stop());
    results_printer.print(warmup, "tajimas_d", "sfkit_dag", trees_file, "tajimas_d", tajimas_d, "1", iteration);

    memory_usage.start();
    timer.start();

    // --- Benchmark computing the FST. ---
    // tskit's implementation is in Python; see Tajima's D above.
    memory_usage.start();
    timer.start();

    double const sfkit_fst = dag_succinct_forest.fst(sample_set_1, sample_set_2);
    do_not_optimize(sfkit_fst);

    log_time(warmup, "fst", "sfkit_dag", timer.stop());
    log_mem(warmup, "fst", "sfkit_dag", memory_usage.stop());
    results_printer.print(warmup, "fst", "sfkit_dag", trees_file, "fst", sfkit_fst, "1", iteration);
}
