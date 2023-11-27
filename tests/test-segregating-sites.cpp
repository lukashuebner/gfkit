
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <tskit.h>

#include "sfkit/SuccinctForest.hpp"
#include "sfkit/assertion_levels.hpp"
#include "sfkit/dag/DAGCompressedForest.hpp"
#include "sfkit/dag/DAGForestCompressor.hpp"
#include "sfkit/stats/AlleleFrequencySpectrum.hpp"
#include "sfkit/tskit/tskit.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch;
using namespace ::Catch::Matchers;

// This test case is taken from the tskit test suite (the only test case for the AFS in there that checks values).
TEST_CASE("Segregating Sites tskit example", "[SegregatingSites]") {
    tsk_treeseq_t tskit_tree_sequence;
    tsk_id_t      samples[]        = {0, 1, 2, 3};
    tsk_size_t    sample_set_sizes = 4;

    double reference_num_seg_sites;
    int    ret;

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

    ret = tsk_treeseq_segregating_sites(
        &tskit_tree_sequence,
        1,
        &sample_set_sizes,
        samples,
        0,
        NULL,
        TSK_STAT_SITE,
        &reference_num_seg_sites
    );
    REQUIRE(ret == 0);
    REQUIRE(reference_num_seg_sites == Approx(3.0).epsilon(1e-6));

    // Test our wrapper around tskit
    sfkit::tskit::TSKitTreeSequence tree_sequence(tskit_tree_sequence); // Takes ownership
    CHECK(tree_sequence.num_segregating_sites() == Approx(reference_num_seg_sites).epsilon(1e-6));

    // Test our implementation on the compressed forest.
    sfkit::SuccinctForest forest(std::move(tree_sequence));

    CHECK(forest.num_segregating_sites() == Approx(reference_num_seg_sites).epsilon(1e-6));

    // Do not free the tskit tree sequence, as we transferred ownershop to  sfkit_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}
