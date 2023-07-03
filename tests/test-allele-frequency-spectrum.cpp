#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <tskit.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/load/compressed-forest-serialization.hpp"
#include "tdt/sequence/allele-frequency-spectrum.hpp"
#include "tdt/tskit.hpp"
#include "tskit-testlib/testlib.hpp"
#include "tdt/graph/compressed-forest.hpp"
#include "tdt/load/forest-compressor.hpp"

using namespace ::Catch::Matchers;

// This test case is taken from the tskit test suite (the only test case for the AFS in there that checks values).
TEST_CASE("AFS example multi tree no back no recurrent", "[AlleleFrequencySpectrum]") {
    tsk_treeseq_t tskit_tree_sequence;
    tsk_id_t      samples[]          = {0, 1, 2, 3};
    tsk_size_t    sample_set_sizes[] = {4, 0};
    double        ref_result[5];
    int           ret;

    tsk_treeseq_from_text(
        &tskit_tree_sequence,
        10,
        paper_ex_nodes,
        paper_ex_edges,
        NULL,
        paper_ex_sites,
        paper_ex_mutations,
        paper_ex_individuals,
        NULL,
        0
    );

    // We have two singletons and one tripleton

    // Unpolarised AFS (not implemented in tdt yet)
    // ret = tsk_treeseq_allele_frequency_spectrum(&ts, 1, sample_set_sizes, samples, 0, NULL, 0, result);
    // CHECK(ret == 0);
    // CHECK(result[0] == 0);
    // CHECK(result[1] == 3.0);
    // CHECK(result[2] == 0);

    // Polarised AFS
    ret = tsk_treeseq_allele_frequency_spectrum(
        &tskit_tree_sequence,
        1,
        sample_set_sizes,
        samples,
        0,
        NULL,
        TSK_STAT_POLARISED,
        ref_result
    );
    CHECK(ret == 0);
    CHECK(ref_result[0] == 0);
    CHECK(ref_result[1] == 2.0);
    CHECK(ref_result[2] == 0);
    CHECK(ref_result[3] == 1.0);
    CHECK(ref_result[4] == 0);

    // Test our method of calling tskit's AFS with our wrapper.
    TSKitTreeSequence tdt_tree_sequence(tskit_tree_sequence); // Takes ownership
    auto              wrapper_result = tdt_tree_sequence.allele_frequency_spectrum();
    CHECK_THAT(wrapper_result, RangeEquals(ref_result));

    ForestCompressor       forest_compressor(tdt_tree_sequence);
    GenomicSequenceFactory sequence_store_factory(tdt_tree_sequence);
    CompressedForest       compressed_forest = forest_compressor.compress(sequence_store_factory);

    auto sequence_store = sequence_store_factory.move_storage();

    AlleleFrequencySpectrum<PerfectNumericHasher> afs(sequence_store, compressed_forest);

    CHECK(afs.num_samples() == 4);
    CHECK(afs.num_sites() == 3);

    CHECK(afs.frequency(0) == 0u);
    CHECK(afs.frequency(1) == 2u);
    CHECK(afs.frequency(2) == 0u);
    CHECK(afs.frequency(3) == 1u);
    CHECK(afs.frequency(4) == 0u);
    // tskit stores 0's in the lowest and highest position of the AFS for unwindowed AFS, we don't.
    CHECK_THAT(
        std::span(afs).subspan(1, afs.num_samples() - 1),
        RangeEquals(std::span(ref_result).subspan(1, afs.num_samples() - 1))
    );

    // Do not free the tskit tree sequence, as we transferred ownershop to  tdt_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}

// This /dataset/ is taken from the tskit test suite -- they don't use it as a unit test for the AFS.
TEST_CASE("AFS example single tree back recurrent", "[AlleleFrequencySpectrum]") {
    tsk_treeseq_t tskit_tree_sequence;
    tsk_id_t      samples[]          = {0, 1, 2, 3};
    tsk_size_t    sample_set_sizes[] = {4, 0};
    double        ref_result[5];
    int           ret;

    tsk_treeseq_from_text(
        &tskit_tree_sequence,
        1,
        single_tree_ex_nodes,
        single_tree_ex_edges,
        NULL,
        single_tree_ex_sites,
        single_tree_ex_mutations,
        NULL,
        NULL,
        0
    );
    REQUIRE(tsk_treeseq_get_num_trees(&tskit_tree_sequence) == 1);

    // Unpolarised AFS not implemented in tdt yet

    // Polarised AFS
    ret = tsk_treeseq_allele_frequency_spectrum(
        &tskit_tree_sequence,
        1,
        sample_set_sizes,
        samples,
        0,
        NULL,
        TSK_STAT_POLARISED,
        ref_result
    );
    CHECK(ret == 0);
    CHECK(ref_result[0] == 0);
    CHECK(ref_result[1] == 2.0);
    CHECK(ref_result[2] == 0);
    CHECK(ref_result[3] == 0);
    CHECK(ref_result[4] == 0); // tskit considers all-sites-derived as no-site-derived

    // Test our method of calling tskit's AFS with our wrapper.
    TSKitTreeSequence tdt_tree_sequence(tskit_tree_sequence); // Takes ownership
    auto              wrapper_result = tdt_tree_sequence.allele_frequency_spectrum();
    CHECK_THAT(wrapper_result, RangeEquals(ref_result));

    ForestCompressor       forest_compressor(tdt_tree_sequence);
    GenomicSequenceFactory sequence_store_factory(tdt_tree_sequence);
    CompressedForest       compressed_forest = forest_compressor.compress(sequence_store_factory);
    GenomicSequence        sequence_store    = sequence_store_factory.move_storage();

    AlleleFrequencySpectrum<PerfectNumericHasher> afs(sequence_store, compressed_forest);
    CHECK(afs.num_samples() == 4);
    CHECK(afs.num_sites() == 3);

    CHECK(afs.frequency(0) == 0u);
    CHECK(afs.frequency(1) == 2u);
    CHECK(afs.frequency(2) == 0u);
    CHECK(afs.frequency(3) == 0u);
    CHECK(afs.frequency(4) == 1u);
    // tskit stores 0's in the lowest and highest position of the AFS for unwindowed AFS, we don't.
    CHECK_THAT(
        std::span(afs).subspan(1, afs.num_samples() - 1),
        RangeEquals(std::span(ref_result).subspan(1, afs.num_samples() - 1))
    );

    // Do not free the tskit tree sequence, as we transferred ownershop to  tdt_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}

TEST_CASE("AFS example multiple derived states", "[AlleleFrequencySpectrum]") {
    tsk_treeseq_t tskit_tree_sequence;
    tsk_id_t      samples[]          = {0, 1, 2, 3};
    tsk_size_t    sample_set_sizes[] = {4, 0};
    double        ref_result[5];
    int           ret;

    tsk_treeseq_from_text(
        &tskit_tree_sequence,
        10,
        multi_derived_states_nodes,
        multi_derived_states_edges,
        NULL,
        multi_derived_states_sites,
        multi_derived_states_mutations,
        multi_derived_states_individuals,
        NULL,
        0
    );

    // We have two singletons and one tripleton

    // Unpolarised AFS (not implemented in tdt yet)
    // ret = tsk_treeseq_allele_frequency_spectrum(&ts, 1, sample_set_sizes, samples, 0, NULL, 0, result);
    // CHECK(ret == 0);
    // CHECK(result[0] == 0);
    // CHECK(result[1] == 3.0);
    // CHECK(result[2] == 0);

    // Polarised AFS
    ret = tsk_treeseq_allele_frequency_spectrum(
        &tskit_tree_sequence,
        1,
        sample_set_sizes,
        samples,
        0,
        NULL,
        TSK_STAT_POLARISED,
        ref_result
    );
    CHECK(ret == 0);
    CHECK(ref_result[0] == 0);
    CHECK(ref_result[1] == 3.0);
    CHECK(ref_result[2] == 1.0);
    CHECK(ref_result[3] == 0);
    CHECK(ref_result[4] == 0);

    // Test our method of calling tskit's AFS with our wrapper.
    TSKitTreeSequence tdt_tree_sequence(tskit_tree_sequence); // Takes ownership
    auto              wrapper_result = tdt_tree_sequence.allele_frequency_spectrum();
    CHECK_THAT(wrapper_result, RangeEquals(ref_result));

    ForestCompressor       forest_compressor(tdt_tree_sequence);
    GenomicSequenceFactory sequence_store_factory(tdt_tree_sequence);
    CompressedForest       compressed_forest = forest_compressor.compress(sequence_store_factory);
    GenomicSequence        sequence_store    = sequence_store_factory.move_storage();

    AlleleFrequencySpectrum<PerfectNumericHasher> afs(sequence_store, compressed_forest);
    CHECK(afs.num_samples() == 4);
    CHECK(afs.num_sites() == 3);

    // TODO Manually check these values!
    CHECK(afs.frequency(0) == 0u);
    CHECK(afs.frequency(1) == 3u);
    CHECK(afs.frequency(2) == 1u);
    CHECK(afs.frequency(3) == 0u);
    CHECK(afs.frequency(4) == 0u);
    // tskit stores 0's in the lowest and highest position of the AFS for unwindowed AFS, we don't.
    CHECK_THAT(
        std::span(afs).subspan(1, afs.num_samples() - 1),
        RangeEquals(std::span(ref_result).subspan(1, afs.num_samples() - 1))
    );

    // Do not free the tskit tree sequence, as we transferred ownershop to  tdt_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}
TEST_CASE("AFS example multi tree back recurrent", "[AlleleFrequencySpectrum]") {
    tsk_treeseq_t tskit_tree_sequence;
    tsk_id_t      samples[]          = {0, 1, 2, 3};
    tsk_size_t    sample_set_sizes[] = {4, 0};
    double        ref_result[5];
    int           ret;

    tsk_treeseq_from_text(
        &tskit_tree_sequence,
        10,
        multi_tree_back_recurrent_nodes,
        multi_tree_back_recurrent_edges,
        NULL,
        multi_tree_back_recurrent_sites,
        multi_tree_back_recurrent_mutations,
        multi_tree_back_recurrent_individuals,
        NULL,
        0
    );

    // We have two singletons and one tripleton

    // Unpolarised AFS (not implemented in tdt yet)
    // ret = tsk_treeseq_allele_frequency_spectrum(&ts, 1, sample_set_sizes, samples, 0, NULL, 0, result);
    // CHECK(ret == 0);
    // CHECK(result[0] == 0);
    // CHECK(result[1] == 3.0);
    // CHECK(result[2] == 0);

    // Polarised AFS
    ret = tsk_treeseq_allele_frequency_spectrum(
        &tskit_tree_sequence,
        1,
        sample_set_sizes,
        samples,
        0,
        NULL,
        TSK_STAT_POLARISED,
        ref_result
    );
    CHECK(ret == 0);
    CHECK(ref_result[0] == 0);
    CHECK(ref_result[1] == 1.0);
    CHECK(ref_result[2] == 2.0);
    CHECK(ref_result[3] == 0);
    CHECK(ref_result[4] == 0);

    // Test our method of calling tskit's AFS with our wrapper.
    TSKitTreeSequence tdt_tree_sequence(tskit_tree_sequence); // Takes ownership
    auto              wrapper_result = tdt_tree_sequence.allele_frequency_spectrum();
    CHECK_THAT(wrapper_result, RangeEquals(ref_result));

    // TODO Write a utility function to construct the compressed tree + sequence store and use it in all the tests.
    ForestCompressor       forest_compressor(tdt_tree_sequence);
    GenomicSequenceFactory sequence_store_factory(tdt_tree_sequence);
    CompressedForest       compressed_forest = forest_compressor.compress(sequence_store_factory);
    GenomicSequence        sequence_store    = sequence_store_factory.move_storage();

    AlleleFrequencySpectrum<PerfectNumericHasher> afs(sequence_store, compressed_forest);
    CHECK(afs.num_samples() == 4);
    CHECK(afs.num_sites() == 3);

    CHECK(afs.frequency(0) == 0u);
    CHECK(afs.frequency(1) == 1u);
    CHECK(afs.frequency(2) == 2u);
    CHECK(afs.frequency(3) == 0u);
    CHECK(afs.frequency(4) == 0u);
    // tskit stores 0's in the lowest and highest position of the AFS for unwindowed AFS, we don't.
    CHECK_THAT(
        std::span(afs).subspan(1, afs.num_samples() - 1),
        RangeEquals(std::span(ref_result).subspan(1, afs.num_samples() - 1))
    );

    // Do not free the tskit tree sequence, as we transferred ownershop to  tdt_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}

TEST_CASE("AFS simulated dataset", "[AlleleFrequencySpectrum]") {
    std::vector<std::string> const ts_files = {
        "data/allele-frequency-spectrum-simple-example-0.trees",
        "data/allele-frequency-spectrum-simple-example-1.trees",
        "data/allele-frequency-spectrum-simple-example-2.trees",
        "data/allele-frequency-spectrum-simple-example-3.trees",
        "data/allele-frequency-spectrum-simple-example-4.trees",
        "data/allele-frequency-spectrum-simple-example-6.trees",
    };
    auto const& ts_file = GENERATE_REF(from_range(ts_files));

    TSKitTreeSequence tree_sequence(ts_file);
    auto              reference_afs = tree_sequence.allele_frequency_spectrum();

    ForestCompressor       forest_compressor(tree_sequence);
    GenomicSequenceFactory sequence_store_factory(tree_sequence);
    CompressedForest       compressed_forest = forest_compressor.compress(sequence_store_factory);
    GenomicSequence        sequence_store    = sequence_store_factory.move_storage();

    AlleleFrequencySpectrum<PerfectDNAHasher> afs(sequence_store, compressed_forest);

    CHECK(afs.num_samples() == compressed_forest.postorder_edges().num_leaves());
    CHECK(afs.num_sites() == sequence_store.num_sites());

    // tskit stores 0's in the lowest and highest position of the AFS for unwindowed AFS, we don't.
    CHECK_THAT(
        std::span(afs).subspan(1, afs.num_samples() - 1),
        RangeEquals(std::span(reference_afs).subspan(1, afs.num_samples() - 1))
    );
}

void convert(std::string const& trees_file, std::string const& forest_file) {
    TSKitTreeSequence tree_sequence(trees_file);

    ForestCompressor       forest_compressor(tree_sequence);
    GenomicSequenceFactory sequence_store_factory(tree_sequence);
    CompressedForest       forest   = forest_compressor.compress(sequence_store_factory);
    GenomicSequence        sequence = sequence_store_factory.move_storage();

    CompressedForestIO::save(forest_file, forest, sequence);
}

TEST_CASE("AFS from .forest", "[AlleleFrequencySpectrum]") {
    struct Dataset {
        std::string forest_file;
        std::string trees_file;
    };

    std::vector<Dataset> const datasets = {
        {"data/allele-frequency-spectrum-simple-example-0.forest",
         "data/allele-frequency-spectrum-simple-example-0.trees"},
        {"data/allele-frequency-spectrum-simple-example-1.forest",
         "data/allele-frequency-spectrum-simple-example-1.trees"},
        {"data/allele-frequency-spectrum-simple-example-2.forest",
         "data/allele-frequency-spectrum-simple-example-2.trees"},
        {"data/allele-frequency-spectrum-simple-example-3.forest",
         "data/allele-frequency-spectrum-simple-example-3.trees"},
        {"data/allele-frequency-spectrum-simple-example-4.forest",
         "data/allele-frequency-spectrum-simple-example-4.trees"},
        {"data/allele-frequency-spectrum-simple-example-6.forest",
         "data/allele-frequency-spectrum-simple-example-6.trees"},
    };

    auto const& dataset     = GENERATE_REF(from_range(datasets));
    auto const& forest_file = dataset.forest_file;
    auto const& trees_file  = dataset.trees_file;

    convert(trees_file, forest_file);

    CompressedForest forest;
    GenomicSequence  sequence;
    CompressedForestIO::load(forest_file, forest, sequence);

    TSKitTreeSequence tree_sequence(trees_file);
    auto              reference_afs = tree_sequence.allele_frequency_spectrum();

    AlleleFrequencySpectrum<PerfectDNAHasher> afs(sequence, forest);

    CHECK(afs.num_samples() == forest.postorder_edges().num_leaves());
    CHECK(afs.num_sites() == sequence.num_sites());

    // tskit stores 0's in the lowest and highest position of the AFS for unwindowed AFS, we don't.
    CHECK_THAT(
        std::span(afs).subspan(1, afs.num_samples() - 1),
        RangeEquals(std::span(reference_afs).subspan(1, afs.num_samples() - 1))
    );
}
