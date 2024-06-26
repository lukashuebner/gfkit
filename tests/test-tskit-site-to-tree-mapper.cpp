
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <tskit.h>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/dag/DAGForestCompressor.hpp"
#include "sfkit/sequence/TSKitSiteToTreeMapper.hpp"
#include "sfkit/tskit/tskit.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch::Matchers;

using sfkit::graph::NodeId;
using sfkit::samples::SampleId;
using sfkit::sequence::SiteId;
using sfkit::sequence::TSKitSiteToTreeMapper;
using sfkit::tskit::TSKitTreeSequence;

// This test case is taken from the tskit test suite (the only test case for the AFS in there that checks values).
TEST_CASE("TSKitSiteToTreeMapper example multi tree no back no recurrent", "[TSKitSiteToTreeMapper]") {
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

    TSKitTreeSequence sfkit_tree_sequence(std::move(tskit_tree_sequence));
    REQUIRE(sfkit_tree_sequence.is_owning());

    TSKitSiteToTreeMapper site2tree(sfkit_tree_sequence);

    CHECK(site2tree(0) == 0);
    CHECK(site2tree(1) == 1);
    CHECK(site2tree(2) == 2);
}

// This /dataset/ is taken from the tskit test suite -- they don't use it as a unit test for the AFS.
TEST_CASE("TSKitSiteToTreeMapper example single tree back recurrent", "[TSKitSiteToTreeMapper]") {
    tsk_treeseq_t tskit_tree_sequence;

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

    TSKitTreeSequence sfkit_tree_sequence(std::move(tskit_tree_sequence));
    REQUIRE(sfkit_tree_sequence.is_owning());

    TSKitSiteToTreeMapper site2tree(sfkit_tree_sequence);

    CHECK(site2tree(0) == 0);
    CHECK(site2tree(1) == 0);
    CHECK(site2tree(2) == 0);
}

TEST_CASE("TSKitSiteToTreeMapper example multi tree back recurrent", "[TSKitSiteToTreeMapper]") {
    tsk_treeseq_t tskit_tree_sequence;

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

    TSKitTreeSequence sfkit_tree_sequence(std::move(tskit_tree_sequence));
    REQUIRE(sfkit_tree_sequence.is_owning());

    TSKitSiteToTreeMapper site2tree(sfkit_tree_sequence);

    CHECK(site2tree(0) == 0);
    CHECK(site2tree(1) == 1);
    CHECK(site2tree(2) == 2);
}
