#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <tskit.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/graph/compressed-forest.hpp"
#include "tdt/load/forest-compressor.hpp"
#include "tdt/succinct-forest.hpp"
#include "tdt/sequence/allele-frequencies.hpp"
#include "tdt/tskit.hpp"
#include "tdt/utils/literals.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch::Matchers;

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

    SuccinctForest sequence_forest(tskit_tree_sequence); // Takes ownership

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

    // Do not free the tskit tree sequence, as we transferred ownershop to  tdt_tree_sequence now.
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

    SuccinctForest<PerfectNumericHasher> sequence_forest(tskit_tree_sequence); // Takes ownership

    // .allele_frequencies() returns the number of samples in the ancestral state.
    auto const& all_samples = sequence_forest.all_samples();
    auto const  freqs       = sequence_forest.allele_frequencies(all_samples);

    auto freqs_it = freqs.begin();
    REQUIRE(freqs_it != freqs.end());

    using MultiallelicFrequency = AlleleFrequencies<PerfectNumericHasher>::MultiallelicFrequency;
    REQUIRE(std::holds_alternative<MultiallelicFrequency>(*freqs_it));
    MultiallelicFrequency const& freq = std::get<MultiallelicFrequency>(*freqs_it);
    CHECK(freq[0] == 1);
    CHECK(freq[1] == 1);
    CHECK(freq[2] == 2);
    CHECK(freq[3] == 0);

    freqs_it++;
    CHECK(freqs_it == freqs.end());

    // Do not free the tskit tree sequence, as we transferred ownershop to  tdt_tree_sequence now.
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

    SuccinctForest<PerfectNumericHasher> sequence_forest(tskit_tree_sequence); // Takes ownership

    // .allele_frequencies() returns the number of samples in the ancestral state.
    // Site 0 is multiallelic with 2 1 1 0 derived states.
    // Site 1 is biallelic with 1 derived state.
    // Site 2 is biallelic with 2 derived states.
    using MultiallelicFrequency = typename AlleleFrequencies<PerfectNumericHasher>::MultiallelicFrequency;
    using BiallelicFrequency    = typename AlleleFrequencies<PerfectNumericHasher>::BiallelicFrequency;
    using AlleleFrequency       = typename AlleleFrequencies<PerfectNumericHasher>::AlleleFrequency;
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

    // Do not free the tskit tree sequence, as we transferred ownershop to  tdt_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}
