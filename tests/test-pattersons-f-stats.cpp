
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

using sfkit::SuccinctForest;
using sfkit::graph::NodeId;
using sfkit::samples::SampleId;
using sfkit::samples::SampleSet;
using sfkit::sequence::PerfectNumericHasher;
using sfkit::sequence::SiteId;
using sfkit::tskit::TSKitTreeSequence;

TEST_CASE("Patterson's f{2,3,4} tskit examples", "[PattersonsFStats]") {
    struct Dataset {
        char const* name;
        char const* nodes;
        char const* edges;
        char const* sites;
        int const   sequence_len;
        char const* mutations;
        char const* individuals;
    };

    std::vector<Dataset> datasets{
        {"paper_ex", paper_ex_nodes, paper_ex_edges, paper_ex_sites, 10, paper_ex_mutations, paper_ex_individuals},
        {"single_tree_ex",
         single_tree_ex_nodes,
         single_tree_ex_edges,
         single_tree_ex_sites,
         1,
         single_tree_ex_mutations,
         NULL},
        {"multi_tree_back_reccurent",
         multi_tree_back_recurrent_nodes,
         multi_tree_back_recurrent_edges,
         multi_tree_back_recurrent_sites,
         10,
         multi_tree_back_recurrent_mutations,
         multi_tree_back_recurrent_individuals},
        {"multi_derived_states",
         multi_derived_states_nodes,
         multi_derived_states_edges,
         multi_derived_states_sites,
         10,
         multi_derived_states_mutations,
         multi_derived_states_individuals}};

    Dataset const& dataset = GENERATE_REF(from_range(datasets));

    tsk_treeseq_t tskit_tree_sequence;

    // f4 setup
    tsk_id_t   f4_samples[]          = {0, 1, 2, 3};
    tsk_size_t f4_sample_set_sizes[] = {1, 1, 1, 1};
    tsk_id_t   f4_set_indexes[]      = {0, 1, 2, 3};
    auto const f4_sample_set_1       = SampleSet(4).add(0);
    auto const f4_sample_set_2       = SampleSet(4).add(1);
    auto const f4_sample_set_3       = SampleSet(4).add(2);
    auto const f4_sample_set_4       = SampleSet(4).add(3);

    // f3 setup
    tsk_id_t   f3_samples[]          = {0, 1, 2, 3};
    tsk_size_t f3_sample_set_sizes[] = {2, 1, 1};
    tsk_id_t   f3_set_indexes[]      = {0, 1, 2};
    auto const f3_sample_set_1       = SampleSet(4).add(0).add(1);
    auto const f3_sample_set_2       = SampleSet(4).add(2);
    auto const f3_sample_set_3       = SampleSet(4).add(3);

    // f2 setup
    tsk_id_t   f2_samples[]          = {0, 1, 2, 3};
    tsk_size_t f2_sample_set_sizes[] = {2, 2};
    tsk_id_t   f2_set_indexes[]      = {0, 1};
    auto const f2_sample_set_1       = SampleSet(4).add(0).add(1);
    auto const f2_sample_set_2       = SampleSet(4).add(2).add(3);

    int ret;
    tsk_treeseq_from_text(
        &tskit_tree_sequence,
        dataset.sequence_len,
        dataset.nodes,
        dataset.edges,
        NULL,
        dataset.sites,
        dataset.mutations,
        dataset.individuals,
        NULL,
        0
    );

    // Test the reference implementation (tskit)
    double reference_f4;
    ret = tsk_treeseq_f4(
        &tskit_tree_sequence,
        4,
        f4_sample_set_sizes,
        f4_samples,
        1,
        f4_set_indexes,
        0,
        NULL,
        TSK_STAT_SITE,
        &reference_f4
    );
    REQUIRE(ret == 0);

    double reference_f3;
    ret = tsk_treeseq_f3(
        &tskit_tree_sequence,
        3,
        f3_sample_set_sizes,
        f3_samples,
        1,
        f3_set_indexes,
        0,
        NULL,
        TSK_STAT_SITE,
        &reference_f3
    );
    REQUIRE(ret == 0);

    double reference_f2;
    ret = tsk_treeseq_f2(
        &tskit_tree_sequence,
        2,
        f2_sample_set_sizes,
        f2_samples,
        1,
        f2_set_indexes,
        0,
        NULL,
        TSK_STAT_SITE,
        &reference_f2
    );
    REQUIRE(ret == 0);

    // Test our wrapper around tskit
    TSKitTreeSequence tree_sequence(tskit_tree_sequence); // Takes ownership of tskit_tree_sequence
    CHECK(
        tree_sequence.f4(f4_sample_set_1, f4_sample_set_2, f4_sample_set_3, f4_sample_set_4)
        == Approx(reference_f4).epsilon(1e-6)
    );
    CHECK(tree_sequence.f3(f3_sample_set_1, f3_sample_set_2, f3_sample_set_3) == Approx(reference_f3).epsilon(1e-6));
    CHECK(tree_sequence.f2(f2_sample_set_1, f2_sample_set_2) == Approx(reference_f2).epsilon(1e-6));

    // Test our implementation on the compressed forest.
    SuccinctForest<PerfectNumericHasher> forest(std::move(tree_sequence));

    CHECK(
        forest.f4(f4_sample_set_1, f4_sample_set_2, f4_sample_set_3, f4_sample_set_4)
        == Approx(reference_f4).epsilon(1e-6)
    );

    CHECK(forest.f3(f3_sample_set_1, f3_sample_set_2, f3_sample_set_3) == Approx(reference_f3).epsilon(1e-6));

    CHECK(forest.f2(f2_sample_set_1, f2_sample_set_2) == Approx(reference_f2).epsilon(1e-6));

    // Do not free the tskit tree sequence, as we transferred ownershop to  sfkit_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}
