#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <tskit.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/graph/compressed-forest.hpp"
#include "tdt/load/forest-compressor.hpp"
#include "tdt/sequence-forest.hpp"
#include "tdt/tskit.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch;
using namespace ::Catch::Matchers;

// This test case is taken from the tskit test suite (the only test case for the AFS in there that checks values).
TEST_CASE("Divergence tskit example", "[Divergence]") {
    tsk_treeseq_t tskit_tree_sequence;
    tsk_id_t      samples[]          = {0, 1, 2, 3};
    tsk_size_t    sample_set_sizes[] = {2, 2};
    tsk_id_t      set_indexes[]      = {0, 1};
    SampleSet     sample_set_1(4);
    SampleSet     sample_set_2(4);
    sample_set_1.add(0);
    sample_set_1.add(1);
    sample_set_2.add(2);
    sample_set_2.add(3);

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

    // Test the reference implementation (tskit)
    double reference_divergence;

    auto ret = tsk_treeseq_divergence(
        &tskit_tree_sequence,
        2,
        sample_set_sizes,
        samples,
        1,
        set_indexes,
        0,
        NULL,
        TSK_STAT_SITE,
        &reference_divergence
    );
    REQUIRE(ret == 0);
    // REQUIRE(reference_pi == Approx(1.5).epsilon(1e-6));

    // Test our wrapper around tskit
    TSKitTreeSequence tree_sequence(tskit_tree_sequence); // Takes ownership of tskit_tree_sequence
    CHECK(tree_sequence.divergence(sample_set_1, sample_set_2) == Approx(reference_divergence).epsilon(1e-6));

    // Test our implementation on the compressed forest.
    SequenceForest sequence_forest(std::move(tree_sequence));

    CHECK(sequence_forest.divergence(sample_set_1, sample_set_2) == Approx(reference_divergence).epsilon(1e-6));

    // Do not free the tskit tree sequence, as we transferred ownershop to  tdt_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}
