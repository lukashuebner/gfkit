#include <filesystem>
#include <fstream>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <sfkit/include-redirects/cereal.hpp>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/load/CompressedForestIO.hpp"
#include "sfkit/load/ForestCompressor.hpp"
#include "sfkit/sequence/AlleleFrequencySpectrum.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"
#include "sfkit/tskit.hpp"

using namespace ::Catch::Matchers;
using Catch::Approx;

std::string const ARCHIVE_FILE_NAME = "tmp-test-f5f515340fa29c848db5ed746253f571c6c791bb.forest";

TEST_CASE("CompressedForest/GenomicSequenceStorage Serialization", "[Serialization]") {
    std::vector<std::string> const ts_files = {
        "data/test-sarafina.trees",
        "data/test-scar.trees",
        "data/test-shenzi.trees",
        "data/test-banzai.trees",
        "data/test-ed.trees",
        "data/test-simba.trees",
    };
    auto const& ts_file = GENERATE_REF(from_range(ts_files));

    TSKitTreeSequence tree_sequence(ts_file);

    ForestCompressor       forest_compressor(tree_sequence);
    GenomicSequenceFactory sequence_factory(tree_sequence);
    CompressedForest       forest   = forest_compressor.compress(sequence_factory);
    GenomicSequence        sequence = sequence_factory.move_storage();

    // Serialize and deserialize the compressed forest and genome sequence storage
    CompressedForestIO::save(ARCHIVE_FILE_NAME, forest, sequence);

    CompressedForest forest_deserialized;
    GenomicSequence  sequence_deserialized;
    CompressedForestIO::load(ARCHIVE_FILE_NAME, forest_deserialized, sequence_deserialized);

    std::filesystem::remove(ARCHIVE_FILE_NAME);

    // CompressedForest
    CHECK(forest_deserialized.num_nodes() == forest.num_nodes());
    CHECK(forest_deserialized.num_samples() == forest.num_samples());
    CHECK(forest_deserialized.num_trees() == forest.num_trees());
    CHECK(forest_deserialized.num_edges() == forest.num_edges());
    CHECK(forest_deserialized.num_roots() == forest.num_roots());
    CHECK(forest_deserialized.roots() == forest.roots());

    auto const all_samples = forest.all_samples();
    auto const num_samples_below_deserialized =
        NumSamplesBelowFactory::build(forest_deserialized.postorder_edges(), all_samples);
    auto const num_samples_below = NumSamplesBelowFactory::build(forest.postorder_edges(), all_samples);
    for (NodeId node_id = 0; node_id < forest.num_nodes(); ++node_id) {
        CHECK(forest_deserialized.is_sample(node_id) == forest.is_sample(node_id));
        CHECK(num_samples_below_deserialized(node_id) == num_samples_below(node_id));
    }

    auto cf_dag              = forest_deserialized.postorder_edges();
    auto cf_dag_deserialized = forest_deserialized.postorder_edges();
    CHECK(cf_dag_deserialized.num_nodes() == cf_dag.num_nodes());
    CHECK(cf_dag_deserialized.num_edges() == cf_dag.num_edges());
    CHECK(cf_dag_deserialized.num_roots() == cf_dag.num_roots());
    CHECK(cf_dag_deserialized.num_trees() == cf_dag.num_trees());
    CHECK(cf_dag_deserialized.roots() == cf_dag.roots());
    CHECK(cf_dag_deserialized.leaves() == cf_dag.leaves());
    CHECK(cf_dag_deserialized.directed() == cf_dag.directed());
    CHECK(cf_dag_deserialized.traversal_order() == cf_dag.traversal_order());
    CHECK(cf_dag_deserialized.is_inorder() == cf_dag.is_inorder());
    CHECK(cf_dag_deserialized.is_postorder() == cf_dag.is_postorder());
    CHECK(cf_dag_deserialized.is_preorder() == cf_dag.is_preorder());
    CHECK(cf_dag_deserialized.is_levelorder() == cf_dag.is_levelorder());
    CHECK(cf_dag_deserialized.is_unordered() == cf_dag.is_unordered());
    CHECK(cf_dag_deserialized.check_postorder() == cf_dag.check_postorder());
    CHECK_THAT(cf_dag_deserialized, RangeEquals(cf_dag));

    // GenomicSequenceStore
    CHECK(sequence_deserialized.num_sites() == sequence.num_sites());
    CHECK(sequence_deserialized.num_mutations() == sequence.num_mutations());

    for (SiteId site_id = 0; site_id < sequence.num_sites(); ++site_id) {
        CHECK(sequence_deserialized.ancestral_state(site_id) == sequence.ancestral_state(site_id));
        CHECK(sequence_deserialized[site_id] == sequence[site_id]);
        CHECK_THAT(sequence_deserialized.mutations_at_site(site_id), RangeEquals(sequence.mutations_at_site(site_id)));
    }

    for (size_t mutation_id = 0; mutation_id < sequence.num_mutations(); ++mutation_id) {
        CHECK(sequence_deserialized.mutation_by_id(mutation_id) == sequence.mutation_by_id(mutation_id));
    }

    CHECK(sequence_deserialized.mutations_are_sorted_by_site() == sequence.mutations_are_sorted_by_site());

    // High-level operations
    SuccinctForest sf(std::move(tree_sequence), std::move(forest), std::move(sequence));

    tree_sequence = TSKitTreeSequence(ts_file);
    SuccinctForest sf_deserialized(
        std::move(tree_sequence),
        std::move(forest_deserialized),
        std::move(sequence_deserialized)
    );

    auto const afs              = sf.allele_frequency_spectrum(sf.all_samples());
    auto const afs_deserialized = sf_deserialized.allele_frequency_spectrum(sf_deserialized.all_samples());

    CHECK(afs.num_samples() == afs_deserialized.num_samples());
    CHECK_THAT(afs, RangeEquals(afs_deserialized));

    // The original and the deserialized forest should have the same sample ids.
    SampleSet sample_set_1(sf.num_samples());
    SampleSet sample_set_2(sf.num_samples());
    bool      flip = false;
    for (SampleId sample: sf.all_samples()) {
        if (flip) {
            sample_set_1.add(sample);
        } else {
            sample_set_2.add(sample);
        }
        flip = !flip;
    }
    CHECK(sf.diversity() == Approx(sf_deserialized.diversity()).epsilon(1e-4));
    CHECK(sf.num_segregating_sites() == sf_deserialized.num_segregating_sites());
    CHECK(
        sf.divergence(sample_set_1, sample_set_2)
        == Approx(sf_deserialized.divergence(sample_set_1, sample_set_2)).epsilon(1e-4)
    );
    CHECK(sf.tajimas_d() == Approx(sf_deserialized.tajimas_d()).epsilon(1e-4));
    CHECK(sf.fst(sample_set_1, sample_set_2) == Approx(sf_deserialized.fst(sample_set_1, sample_set_2)).epsilon(1e-4));
}

void convert(std::string const& trees_file, std::string const& forest_file) {
    TSKitTreeSequence      tree_sequence(trees_file);
    ForestCompressor       forest_compressor(tree_sequence);
    GenomicSequenceFactory sequence_store_factory(tree_sequence);
    CompressedForest       forest   = forest_compressor.compress(sequence_store_factory);
    GenomicSequence        sequence = sequence_store_factory.move_storage();

    CompressedForestIO::save(forest_file, forest, sequence);
}

TEST_CASE("Statistics on .forest files", "[Serialization]") {
    struct Dataset {
        std::string forest_file;
        std::string trees_file;
    };

    // TODO Rename these sample datasets
    // TODO Generate more sample datasets
    std::vector<Dataset> const datasets = {
        {"data/test-sarafina.forest", "data/test-sarafina.trees"},
        {"data/test-scar.forest", "data/test-scar.trees"},
        {"data/test-shenzi.forest", "data/test-shenzi.trees"},
        {"data/test-banzai.forest", "data/test-banzai.trees"},
        {"data/test-ed.forest", "data/test-ed.trees"},
        {"data/test-simba.forest", "data/test-simba.trees"},
    };

    auto const& dataset     = GENERATE_REF(from_range(datasets));
    auto const& forest_file = dataset.forest_file;
    auto const& trees_file  = dataset.trees_file;

    convert(trees_file, forest_file);

    CompressedForest  forest;
    GenomicSequence   sequence;
    TSKitTreeSequence tree_sequence(trees_file);
    CompressedForestIO::load(forest_file, forest, sequence);
    CHECK(forest.num_nodes() > 0);
    CHECK(forest.num_nodes_is_set());

    // --- AFS ---
    auto tskit_afs = tree_sequence.allele_frequency_spectrum();

    AlleleFrequencies<PerfectDNAHasher> const allele_frequencies(forest, sequence, forest.all_samples());
    AlleleFrequencySpectrum const             sfkit_afs(allele_frequencies);

    CHECK(sfkit_afs.num_samples() == forest.postorder_edges().num_leaves());

    // tskit stores 0's in the lowest and highest position of the AFS for unwindowed AFS, we don't.
    for (size_t i = 1; i < sfkit_afs.num_samples() - 1; ++i) {
        if (sfkit_afs[i] != tskit_afs[i]) {
            std::cerr << "ERROR !! AFS mismatch between tskit and sfkit" << std::endl;
            std::cerr << "    " << tskit_afs[i] << " vs. " << sfkit_afs[i] << " (tskit vs. sfkit)" << std::endl;
            std::exit(1);
        }
    }

    // --- Divergence ---
    SuccinctForest sequence_forest(std::move(tree_sequence), std::move(forest), std::move(sequence));

    {
        SampleSet sample_set_1(sequence_forest.num_samples());
        SampleSet sample_set_2(sequence_forest.num_samples());
        bool      flip = false;
        for (SampleId sample: forest.leaves()) {
            if (flip) {
                sample_set_1.add(sample);
            } else {
                sample_set_2.add(sample);
            }
            flip = !flip;
        }

        double const sfkit_divergence = sequence_forest.divergence(sample_set_1, sample_set_2);
        double const tskit_divergence = sequence_forest.tree_sequence().divergence(sample_set_1, sample_set_2);
        CHECK(sfkit_divergence == Approx(tskit_divergence).epsilon(1e-4));
    }

    // --- Diversity ---
    {
        double const sfkit_diversity = sequence_forest.diversity();
        double const tskit_diversity = sequence_forest.tree_sequence().diversity();
        CHECK(sfkit_diversity == Approx(tskit_diversity).epsilon(1e-4));
    }

    // --- f{2,3,4} ---
    constexpr size_t       num_sample_sets = 4;
    std::vector<SampleSet> sample_sets(num_sample_sets, sequence_forest.num_samples());
    size_t                 idx = 0;
    for (SampleId sample: sequence_forest.forest().leaves()) {
        sample_sets[idx].add(sample);
        idx = (idx + 1ul) % num_sample_sets;
    }
    double const sfkit_f4 = sequence_forest.f4(sample_sets[0], sample_sets[1], sample_sets[2], sample_sets[3]);
    double const tskit_f4 =
        sequence_forest.tree_sequence().f4(sample_sets[0], sample_sets[1], sample_sets[2], sample_sets[3]);
    CHECK(sfkit_f4 == Approx(tskit_f4).epsilon(1e-4));

    double const sfkit_f3 = sequence_forest.f3(sample_sets[0], sample_sets[1], sample_sets[2]);
    double const tskit_f3 = sequence_forest.tree_sequence().f3(sample_sets[0], sample_sets[1], sample_sets[2]);
    CHECK(sfkit_f3 == Approx(tskit_f3).epsilon(1e-4));

    double const sfkit_f2 = sequence_forest.f2(sample_sets[0], sample_sets[1]);
    double const tskit_f2 = sequence_forest.tree_sequence().f2(sample_sets[0], sample_sets[1]);
    CHECK(sfkit_f2 == Approx(tskit_f2).epsilon(1e-4));

    // Benchmark computing the number of segregating sites
    SiteId const sfkit_num_seg_sites = sequence_forest.num_segregating_sites();
    double const tskit_num_seg_sites = sequence_forest.tree_sequence().num_segregating_sites();
    CHECK(sfkit_num_seg_sites == Approx(tskit_num_seg_sites).epsilon(1e-4));

    // Check some properties of the compressed forest
    auto const dag = forest.postorder_edges();
    CHECK(dag.check_postorder());
    CHECK(dag.check_no_duplicate_edges());
}
