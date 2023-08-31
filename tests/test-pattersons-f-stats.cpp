
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
#include "tdt/sequence/allele-frequency-spectrum.hpp"
#include "tdt/succinct-forest.hpp"
#include "tdt/tskit.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch;
using namespace ::Catch::Matchers;

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
    tsk_id_t      samples[]          = {0, 1, 2, 3};
    tsk_size_t    sample_set_sizes[] = {1, 1, 1, 1};
    tsk_id_t      set_indexes[]      = {0, 1, 2, 3};
    SampleSet     sample_set_1(4);
    SampleSet     sample_set_2(4);
    SampleSet     sample_set_3(4);
    SampleSet     sample_set_4(4);
    sample_set_1.add(0);
    sample_set_2.add(1);
    sample_set_3.add(2);
    sample_set_4.add(3);

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
        sample_set_sizes,
        samples,
        1,
        set_indexes,
        0,
        NULL,
        TSK_STAT_SITE,
        &reference_f4
    );
    REQUIRE(ret == 0);

    // Test our wrapper around tskit
    TSKitTreeSequence tree_sequence(tskit_tree_sequence); // Takes ownership of tskit_tree_sequence
    CHECK(
        tree_sequence.f4(sample_set_1, sample_set_2, sample_set_3, sample_set_4) == Approx(reference_f4).epsilon(1e-6)
    );

    // Test our implementation on the compressed forest.
    SuccinctForest<PerfectNumericHasher> sequence_forest(std::move(tree_sequence));

    CHECK(
        sequence_forest.f4(sample_set_1, sample_set_2, sample_set_3, sample_set_4) == Approx(reference_f4).epsilon(1e-6)
    );

    // Do not free the tskit tree sequence, as we transferred ownershop to  tdt_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}
