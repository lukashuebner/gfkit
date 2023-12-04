#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <tskit.h>

#include "sfkit/SuccinctForest.hpp"
#include "sfkit/assertion_levels.hpp"
#include "sfkit/dag/DAGCompressedForest.hpp"
#include "sfkit/dag/DAGForestCompressor.hpp"
#include "sfkit/sequence/AlleleFrequencies.hpp"
#include "sfkit/tskit/tskit.hpp"
#include "sfkit/utils/literals.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch::Matchers;
using namespace sfkit;
using sfkit::utils::operator""_uc;

// This test case is taken from the tskit test suite (the only test case for the AFS in there that checks values).
TEST_CASE("AlleleFrequencies biallelic example", "[AlleleFrequencies]") {
    tsk_treeseq_t tskit_tree_sequence;

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

    DAGSuccinctForestNumeric sequence_forest(tskit_tree_sequence);

    // .allele_frequencies() returns the number of samples in the ancestral state.
    std::vector<uint64_t> expected    = {3, 3, 1};
    auto const&           all_samples = sequence_forest.all_samples();
    auto const            freqs       = sequence_forest.allele_frequencies(all_samples);

    size_t idx = 0;
    freqs.visit(
        [&expected, &idx](auto&& state) {
            CHECK(state.num_ancestral() == expected[idx]);
            idx++;
        },
        []([[maybe_unused]] auto&& state) {
            // No multiallelic state
            FAIL();
        }
    );

    // Do not free the tskit tree sequence, as we transferred ownershop to  sfkit_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}

TEST_CASE("AlleleFrequencies multiallelic example", "[AlleleFrequencies]") {
    tsk_treeseq_t tskit_tree_sequence;

    tsk_treeseq_from_text(
        &tskit_tree_sequence,
        1,
        single_tree_multi_derived_states_nodes,
        single_tree_multi_derived_states_edges,
        NULL,
        single_tree_multi_derived_states_sites,
        single_tree_multi_derived_states_mutations,
        NULL,
        NULL,
        0
    );

    DAGSuccinctForestNumeric sequence_forest(tskit_tree_sequence); // Takes ownership

    // .allele_frequencies() returns the number of samples in the ancestral state.
    auto const& all_samples = sequence_forest.all_samples();
    auto const  freqs       = sequence_forest.allele_frequencies(all_samples);

    auto freqs_it = freqs.begin();
    REQUIRE(freqs_it != freqs.end());

    using MultiallelicFrequency = MultiallelicFrequency<PerfectNumericHasher>;
    REQUIRE(std::holds_alternative<MultiallelicFrequency>(*freqs_it));
    MultiallelicFrequency const& freq = std::get<MultiallelicFrequency>(*freqs_it);
    CHECK(freq[0] == 1);
    CHECK(freq[1] == 1);
    CHECK(freq[2] == 2);
    CHECK(freq[3] == 0);

    freqs_it++;
    CHECK(freqs_it == freqs.end());

    // Do not free the tskit tree sequence, as we transferred ownershop to  sfkit_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}

TEST_CASE("AlleleFrequencies mixed example", "[AlleleFrequencies]") {
    tsk_treeseq_t tskit_tree_sequence;

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

    SuccinctForest<DAGCompressedForest, PerfectNumericHasher> sequence_forest(tskit_tree_sequence); // Takes ownership

    // .allele_frequencies() returns the number of samples in the ancestral state.
    // Site 0 is multiallelic with 2 1 1 0 derived states.
    // Site 1 is biallelic with 1 derived state.
    // Site 2 is biallelic with 2 derived states.
    using MultiallelicFrequency = MultiallelicFrequency<PerfectNumericHasher>;
    using AlleleFrequency       = AlleleFrequency<PerfectNumericHasher>;
    std::vector<AlleleFrequency> expected(std::initializer_list<AlleleFrequency>{
        MultiallelicFrequency(static_cast<unsigned char>('0'), 2_uc, 1_uc, 1_uc, 0_uc),
        BiallelicFrequency{3_uc},
        BiallelicFrequency{2_uc}});

    auto const& all_samples = sequence_forest.all_samples();
    auto const  freqs       = sequence_forest.allele_frequencies(all_samples);

    size_t idx = 0;
    freqs.visit(
        [&expected, &idx](auto&& state) {
            CHECK(state == std::get<BiallelicFrequency>(expected[idx]));
            idx++;
        },
        [&expected, &idx](auto&& state) {
            // No multiallelic state
            CHECK(state == std::get<MultiallelicFrequency>(expected[idx]));
            idx++;
        }
    );

    // Do not free the tskit tree sequence, as we transferred ownershop to  sfkit_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}

TEST_CASE("AlleleFrequencies single-sample sample sets", "[AlleleFrequencies]") {
    tsk_treeseq_t tskit_tree_sequence;

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

    SuccinctForest<DAGCompressedForest, PerfectNumericHasher> forest(tskit_tree_sequence); // Takes ownership

    using MultiallelicFrequency = MultiallelicFrequency<PerfectNumericHasher>;
    using AlleleFrequency       = AlleleFrequency<PerfectNumericHasher>;

    auto const sample_0 = SampleSet(4).add(0);
    auto const sample_1 = SampleSet(4).add(1);
    auto const sample_2 = SampleSet(4).add(2);
    auto const sample_3 = SampleSet(4).add(3);

    auto const freqs_0 = forest.allele_frequencies(sample_0);
    auto const freqs_1 = forest.allele_frequencies(sample_1);
    auto const freqs_2 = forest.allele_frequencies(sample_2);
    auto const freqs_3 = forest.allele_frequencies(sample_3);

    std::vector<AlleleFrequency> expected_0(std::initializer_list<AlleleFrequency>{
        BiallelicFrequency{0_uc},
        BiallelicFrequency{1_uc},
        BiallelicFrequency{1_uc}});
    std::vector<AlleleFrequency> expected_1(std::initializer_list<AlleleFrequency>{
        BiallelicFrequency{1_uc},
        BiallelicFrequency{0_uc},
        BiallelicFrequency{1_uc}});
    std::vector<AlleleFrequency> expected_2(std::initializer_list<AlleleFrequency>{
        BiallelicFrequency{1_uc},
        BiallelicFrequency{1_uc},
        BiallelicFrequency{0_uc}});
    std::vector<AlleleFrequency> expected_3(std::initializer_list<AlleleFrequency>{
        MultiallelicFrequency(static_cast<unsigned char>('0'), 0_uc, 0_uc, 1_uc, 0_uc),
        BiallelicFrequency{1_uc},
        BiallelicFrequency{0_uc}});

    size_t idx = 0;
    freqs_0.visit(
        [&expected_0, &idx](auto&& state) {
            CHECK(state == std::get<BiallelicFrequency>(expected_0[idx]));
            idx++;
        },
        []([[maybe_unused]] auto&& state) { FAIL(); }
    );

    idx = 0;
    freqs_1.visit(
        [&expected_1, &idx](auto&& state) {
            CHECK(state == std::get<BiallelicFrequency>(expected_1[idx]));
            idx++;
        },
        []([[maybe_unused]] auto&& state) { FAIL(); }
    );

    idx = 0;
    freqs_2.visit(
        [&expected_2, &idx](auto&& state) {
            CHECK(state == std::get<BiallelicFrequency>(expected_2[idx]));
            idx++;
        },
        []([[maybe_unused]] auto&& state) { FAIL(); }
    );

    idx = 0;
    freqs_3.visit(
        [&expected_3, &idx](auto&& state) {
            CHECK(state == std::get<BiallelicFrequency>(expected_3[idx]));
            idx++;
        },
        [&expected_3, &idx]([[maybe_unused]] auto&& state) {
            CHECK(state == std::get<MultiallelicFrequency>(expected_3[idx]));
            idx++;
        }
    );

    // Do not free the tskit tree sequence, as we transferred ownershop to  sfkit_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}
