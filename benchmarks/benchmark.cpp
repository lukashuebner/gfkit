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
#include "sfkit/SuccinctForest.hpp"
#include "sfkit/bp/BPForestCompressor.hpp"
#include "sfkit/dag/DAGCompressedForest.hpp"
#include "sfkit/dag/DAGForestCompressor.hpp"
#include "sfkit/io/CompressedForestIO.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"
#include "sfkit/stats/AlleleFrequencySpectrum.hpp"
#include "sfkit/tskit/tskit.hpp"
#include "timer.hpp"

std::vector<sfkit::samples::SampleSet>
distribute_samples_to_sample_sets(sfkit::SampleId const num_samples, sfkit::SampleSetId const num_sample_sets) {
    std::vector<sfkit::samples::SampleSet> sample_sets(num_sample_sets, num_samples);

    size_t idx = 0;
    for (sfkit::SampleId sample = 0; sample < num_samples; sample++) {
        sample_sets[idx].add(sample);
        idx = (idx + 1ul) % num_sample_sets;
    }

    return sample_sets;
}

std::vector<sfkit::SampleId> pick_every_nth_sample(sfkit::SampleId const num_overall_samples, sfkit::SampleId const n) {
    std::vector<sfkit::SampleId> output;
    output.reserve(num_overall_samples / n);

    for (sfkit::SampleId i = 0; i < num_overall_samples; i += n) {
        output.push_back(i);
    }

    return output;
}

std::unordered_set<sfkit::SampleId>
pick_n_samples_at_random(sfkit::SampleId const num_overall_samples, sfkit::SampleId const n) {
    std::random_device                             rand_dev;
    std::mt19937                                   generator(rand_dev());
    std::uniform_int_distribution<sfkit::SampleId> pick_sample_at_random(0, num_overall_samples - 1);

    std::unordered_set<sfkit::SampleId> output;
    output.reserve(n);

    while (output.size() < n) {
        output.insert(pick_sample_at_random(generator));
    }

    return output;
}

std::vector<std::pair<sfkit::SampleId, sfkit::SampleId>>
pick_sample_pairs_uniformly_at_random(sfkit::SampleId const num_samples, sfkit::SampleId const n) {
    std::vector<std::pair<sfkit::SampleId, sfkit::SampleId>> pairs;
    pairs.reserve(n);

    for (sfkit::SampleId i = 0; i < n; i++) {
        auto const samples = pick_n_samples_at_random(num_samples, 2);
        KASSERT(samples.size() == 2, "Expected to pick exactly two samples.", sfkit::assert::light);
        auto const first  = *samples.begin();
        auto const second = *samples.begin().operator++();
        KASSERT(first != second, "Expected to pick two different samples.", sfkit::assert::light);
        pairs.emplace_back(first, second);
    }

    return pairs;
}

void validate_stat(
    float sfkit_value, float tskit_value, std::string const& stat_name, std::string const& variant_name
) {
    static constexpr double FLOAT_EQ_EPS = 1e-4;

    if (sfkit_value != Catch::Approx(tskit_value).epsilon(FLOAT_EQ_EPS)) {
        std::cerr << "ERROR !! " << stat_name << " mismatch between tskit and sfkit (" << variant_name << ")"
                  << std::endl;
        std::cerr << "    " << tskit_value << " vs. " << sfkit_value << " (tskit vs. sfkit)" << std::endl;
        // std::exit(EXIT_FAILURE);
    }
}

void benchmark(
    bool const         warmup,
    uint8_t const      iteration,
    std::string const& trees_file,
    std::string const& dag_forest_file,
    std::string const& bp_forest_file,
    ResultsPrinter&    results_printer
) {
    using SetOfSampleSets = sfkit::sequence::SetOfSampleSets<1>;

    KASSERT(
        dag_forest_file != "" || bp_forest_file != "",
        "Specify either DAG or BP forest files or both.",
        sfkit::assert::light
    );

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
    BenchmarkRunner bench(memory_usage, timer, log_mem, log_time);

    // ------ tskit tree sequence ------
    // --- Benchmark tree sequence loading ---
    bench.start();
    sfkit::tskit::TSKitTreeSequence tree_sequence(trees_file);
    do_not_optimize(tree_sequence);
    bench.stop("load_trees_file", "tskit");

    // --- Build sample sets ---
    auto const two_sample_sets   = distribute_samples_to_sample_sets(tree_sequence.num_samples(), 2);
    auto const three_sample_sets = distribute_samples_to_sample_sets(tree_sequence.num_samples(), 3);
    auto const four_sample_sets  = distribute_samples_to_sample_sets(tree_sequence.num_samples(), 4);

    // --- Benchmark computing the LCA/MRCA ---
    uint16_t const num_lca_queries = 1;
    auto const lca_sample_pairs = pick_sample_pairs_uniformly_at_random(tree_sequence.num_samples(), num_lca_queries);

    bench.start();
    for (auto const& [u, v]: lca_sample_pairs) {
        auto const tskit_lca = tree_sequence.lca(u, v);
        do_not_optimize(tskit_lca);
    }
    bench.stop("lca_pairwise", "tskit");
    // The node IDs of sfkit and tskit are not trivially comparable -> Full verification only in the unit tests.

    // tskit's C interface only supports pair-wise queries. We use the Python interface for querying mor than two
    // samples at once.

    // --- Benchmark computing the AFS ---
    bench.start();
    auto const tskit_afs = tree_sequence.allele_frequency_spectrum();
    do_not_optimize(tskit_afs);
    bench.stop("afs", "tskit");
    auto validate_afs = [&tskit_afs](auto const& afs, std::string const& label) {
        for (size_t i = 1; i < afs.num_samples() - 1; ++i) {
            if (afs[i] != tskit_afs[i]) {
                std::cerr << "ERROR !! AFS mismatch between tskit and " << label << std::endl;
                std::cerr << "    " << tskit_afs[i] << " vs. " << afs[i] << " (tskit vs. " << label << ")" << std::endl;
                // std::exit(EXIT_FAILURE);
            }
        }
    };

    // --- Benchmark computing the divergence ---
    bench.start();
    double const tskit_divergence = tree_sequence.divergence(two_sample_sets[0], two_sample_sets[1]);
    do_not_optimize(tskit_divergence);
    bench.stop("divergence", "tskit");
    auto const validate_divergence = [&tskit_divergence](auto const& divergence, std::string const& label) {
        validate_stat(divergence, tskit_divergence, "divergence", label);
    };

    // --- Benchmark computing Patterson's f{2,3,4} ---
    bench.start();
    double const tskit_f4 =
        tree_sequence.f4(four_sample_sets[0], four_sample_sets[1], four_sample_sets[2], four_sample_sets[3]);
    do_not_optimize(tskit_f4);
    bench.stop("f4", "tskit");
    auto const validate_f4 = [&tskit_f4](auto const& f4, std::string const& label) {
        validate_stat(f4, tskit_f4, "f4", label);
    };

    bench.start();
    double const tskit_f3 = tree_sequence.f3(three_sample_sets[0], three_sample_sets[1], three_sample_sets[2]);
    do_not_optimize(tskit_f3);
    bench.stop("f3", "tskit");
    auto const validate_f3 = [&tskit_f3](auto const& f3, std::string const& label) {
        validate_stat(f3, tskit_f3, "f3", label);
    };

    bench.start();
    double const tskit_f2 = tree_sequence.f2(two_sample_sets[0], two_sample_sets[1]);
    do_not_optimize(tskit_f2);
    bench.stop("f2", "tskit");
    auto const validate_f2 = [&tskit_f2](auto const& f2, std::string const& label) {
        validate_stat(f2, tskit_f2, "f2", label);
    };

    // --- Benchmark computing the diversity ---
    bench.start();
    double const tskit_diversity = tree_sequence.diversity();
    do_not_optimize(tskit_diversity);
    bench.stop("diversity", "tskit");
    auto const validate_diversity = [&tskit_diversity](auto const& diversity, std::string const& label) {
        validate_stat(diversity, tskit_diversity, "diversity", label);
    };

    // --- Benchmark computing the number of segregating sites ---
    bench.start();
    double const tskit_num_segregating_sites = tree_sequence.num_segregating_sites();
    do_not_optimize(tskit_num_segregating_sites);
    bench.stop("num_segregating_sites", "tskit");
    auto const validate_num_segregating_sites =
        [&tskit_num_segregating_sites](auto const& num_segregating_sites, std::string const& label) {
            validate_stat(num_segregating_sites, tskit_num_segregating_sites, "num_segregating_sites", label);
        };

    sfkit::dag::DAGCompressedForest         dag_forest;
    sfkit::sequence::GenomicSequenceFactory dag_sequence_factory(tree_sequence);
    sfkit::sequence::GenomicSequence        dag_sequence;
    std::optional<sfkit::DAGSuccinctForest> dag_succinct_forest;

    sfkit::bp::BPCompressedForest           bp_forest;
    sfkit::sequence::GenomicSequenceFactory bp_sequence_factory(tree_sequence);
    sfkit::sequence::GenomicSequence        bp_sequence;
    std::optional<sfkit::DAGSuccinctForest> bp_succinct_forest;

    // --- Benchmark loading the tree sequence ---
    if (dag_forest_file != "") {
        bench.start();
        sfkit::io::CompressedForestIO::load(dag_forest_file, dag_forest, dag_sequence);
        do_not_optimize(dag_forest);
        do_not_optimize(dag_sequence);
        bench.stop("load_forest_file", "sfkit_dag");
        dag_succinct_forest.emplace(std::move(dag_forest), std::move(dag_sequence));
    }

    if (bp_forest_file != "") {
        bench.start();
        sfkit::io::CompressedForestIO::load(bp_forest_file, bp_forest, bp_sequence);
        do_not_optimize(bp_forest);
        do_not_optimize(bp_sequence);
        bench.stop("load_forest_file", "sfkit_bp");
        bp_succinct_forest.emplace(std::move(dag_forest), std::move(dag_sequence));
    }

    // --- Benchmark the sfkit variants ---
    auto const bench_variant = [&](auto succinct_forest, std::string const& variant_name) {
        // --- Benchmark computing the AFS ---
        bench.start();
        auto const sfkit_afs = succinct_forest.allele_frequency_spectrum();
        do_not_optimize(sfkit_afs);
        bench.stop("afs", variant_name);
        validate_afs(sfkit_afs, variant_name);

        // --- Benchmark computing the LCA/MRCA ---
        if constexpr (std::is_same_v<decltype(succinct_forest), sfkit::DAGSuccinctForest>) {
            bench.start();
            for (auto [u, v]: lca_sample_pairs) {
                auto const sfkit_lca = succinct_forest.lca(u, v);
                do_not_optimize(sfkit_lca);
            }
            bench.stop("lca", variant_name);

            for (uint8_t nth = 2; nth <= 10; ++nth) {
                auto const samples_vec = pick_every_nth_sample(tree_sequence.num_samples(), nth);

                sfkit::SampleSet samples(tree_sequence.num_samples());
                for (auto const sample: samples_vec) {
                    samples.add(sample);
                }

                bench.start();
                auto const sfkit_lca = succinct_forest.lca(samples);
                do_not_optimize(sfkit_lca);
                bench.stop(std::string("lca_") + std::to_string(nth) + "th", "sfkit");
            }
        }

        // --- Benchmark computing the divergence ---
        bench.start();
        double const sfkit_divergence = succinct_forest.divergence(two_sample_sets[0], two_sample_sets[1]);
        do_not_optimize(sfkit_divergence);
        bench.stop("divergence", variant_name);
        validate_divergence(sfkit_divergence, variant_name);

        // --- Benchmark computing Patterson's f{2,3,4} ---
        bench.start();
        double const sfkit_f4 =
            succinct_forest.f4(four_sample_sets[0], four_sample_sets[1], four_sample_sets[2], four_sample_sets[3]);
        do_not_optimize(sfkit_f4);
        bench.stop("f4", variant_name);
        validate_f4(sfkit_f4, variant_name);

        bench.start();
        double const sfkit_f3 = succinct_forest.f3(three_sample_sets[0], three_sample_sets[1], three_sample_sets[2]);
        do_not_optimize(sfkit_f3);
        bench.stop("f3", variant_name);
        validate_f3(sfkit_f3, variant_name);

        bench.start();
        double const sfkit_f2 = succinct_forest.f2(two_sample_sets[0], two_sample_sets[1]);
        do_not_optimize(sfkit_f2);
        bench.stop("f2", variant_name);
        validate_f2(sfkit_f2, variant_name);

        // --- Benchmark computing the diversity ---
        bench.start();
        double const sfkit_diversity = succinct_forest.diversity();
        do_not_optimize(sfkit_diversity);
        bench.stop("diversity", variant_name);
        validate_diversity(sfkit_diversity, variant_name);

        // --- Benchmark computing the number of segregating sites ---
        bench.start();
        double const sfkit_num_segregating_sites = succinct_forest.num_segregating_sites();
        do_not_optimize(sfkit_num_segregating_sites);
        bench.stop("num_segregating_sites", variant_name);
        validate_num_segregating_sites(sfkit_num_segregating_sites, variant_name);

        // --- Benchmark computing Tajima's D. ---
        // tskit does implement this only in Python, not in C. We're using
        // experiments/scripts/benchmark-tskits-tajimas-d.py to measure tskit's performance.
        bench.start();
        double const sfkit_tajimas_d = succinct_forest.tajimas_d();
        do_not_optimize(sfkit_tajimas_d);
        bench.stop("tajimas_d", variant_name);

        results_printer
            .print(warmup, "tajimas_d", variant_name, trees_file, "tajimas_d", sfkit_tajimas_d, "1", iteration);

        // --- Benchmark computing the FST. ---
        // tskit's implementation is in Python; see Tajima's D above.
        bench.start();
        double const sfkit_fst = succinct_forest.fst(two_sample_sets[0], two_sample_sets[1]);
        do_not_optimize(sfkit_fst);
        bench.stop("fst", variant_name);

        results_printer.print(warmup, "fst", variant_name, trees_file, "fst", sfkit_fst, "1", iteration);
    };

    if (dag_succinct_forest) {
        bench_variant(*dag_succinct_forest, "sfkit_dag");
    }
    if (bp_succinct_forest) {
        bench_variant(*bp_succinct_forest, "sfkit_bp");
    }
}
