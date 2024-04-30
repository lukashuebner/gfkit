
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>

#include "mocks/TsToSfMappingExtractor.hpp"
#include "sfkit/SuccinctForest.hpp"
#include "sfkit/assertion_levels.hpp"
#include "sfkit/dag/DAGForestCompressor.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/graph/primitives.hpp"
#include "sfkit/stats/LCA.hpp"
#include "sfkit/utils/concepts.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace Catch::Matchers;

using sfkit::dag::DAGCompressedForest;
using sfkit::dag::DAGForestCompressor;
using sfkit::graph::NodeId;
using sfkit::samples::SampleId;
using sfkit::stats::DAGLowestCommonAncestor;
using sfkit::tskit::TSKitTree;
using sfkit::tskit::TSKitTreeSequence;

std::pair<tsk_id_t, tsk_id_t> pick_2_distinct_samples_at_random(tsk_id_t const num_overall_samples) {
    std::random_device                      rand_dev;
    std::mt19937                            generator(rand_dev());
    std::uniform_int_distribution<tsk_id_t> pick_sample_at_random(0, num_overall_samples - 1);

    tsk_id_t const first = pick_sample_at_random(generator);
    tsk_id_t       second;
    do {
        second = pick_sample_at_random(generator);
    } while (first == second);

    return {first, second};
}

TEST_CASE("CompressedForest::lca()", "[CompressedForest]") {
    //                                      38
    //                                ┏━━━━━━┻━━━━━━━┓
    //                               37              ┃
    //                         ┏━━━━━━┻━━━━━━┓       ┃
    //                        36             ┃       ┃
    //                 ┏━━━━━━━┻━━━━━━━━┓    ┃       ┃
    //                 ┃                ┃   35       ┃
    //                 ┃                ┃   ┏┻┓      ┃
    //                 ┃                ┃   ┃ ┃     34
    //                 ┃                ┃   ┃ ┃    ┏━┻━┓
    //                33                ┃   ┃ ┃    ┃   ┃
    //        ┏━━━━━━━━┻━━━━━━━━┓       ┃   ┃ ┃    ┃   ┃
    //       32                 ┃       ┃   ┃ ┃    ┃   ┃
    //    ┏━━━┻━━━┓             ┃       ┃   ┃ ┃    ┃   ┃
    //   31       ┃             ┃       ┃   ┃ ┃    ┃   ┃
    //  ┏━┻━┓     ┃             ┃       ┃   ┃ ┃    ┃   ┃
    //  ┃   ┃     ┃             ┃       ┃   ┃ ┃   30   ┃
    //  ┃   ┃     ┃             ┃       ┃   ┃ ┃  ┏━┻┓  ┃
    //  ┃   ┃    29             ┃       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┏━━━┻━━━┓         ┃       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃      28         ┃       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┏━━┻━━┓      ┃       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┃     ┃     27       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┃     ┃   ┏━━┻━━┓    ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┃     ┃  26     ┃    ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┃     ┃ ┏━┻━┓   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃   25     ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃  ┏━┻━━┓  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // 24   ┃ ┃  ┃    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // ┏┻┓  ┃ ┃  ┃    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ 23    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┏┻━┓  ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃   ┃   ┃   22   ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃   ┃   ┃  ┏━┻━┓ ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃  21   ┃  ┃   ┃ ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┏━┻┓  ┃  ┃   ┃ ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┃  ┃  ┃ 20   ┃ ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┃  ┃  ┃ ┏┻┓  ┃ ┃ ┃  ┃  ┃  ┃
    // 0 6 10 8 9 11 12 13 4 7 14 19 1 2 18 3 5 15 16 17

    std::vector<std::string> const ts_files = {
        "data/test-sarafina.trees",
        "data/test-timon.trees",
        "data/test-scar.trees",
        "data/test-shenzi.trees",
        "data/test-banzai.trees",
        "data/test-ed.trees",
        "data/test-simba.trees",
    };
    auto const&         ts_file = GENERATE_REF(from_range(ts_files));
    TSKitTreeSequence   tree_sequence(ts_file);
    DAGForestCompressor forest_compressor(tree_sequence);

    REQUIRE(tree_sequence.num_trees() >= 1);

    Ts2SfMappingExtractor ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    DAGCompressedForest   forest = forest_compressor.compress(ts_2_sf_node);
    CHECK(ts_2_sf_node.finalize_called());
    CHECK(ts_2_sf_node.process_mutations_callcnt() == tree_sequence.num_trees());

    constexpr uint32_t n_trials = 100;
    for (uint32_t trial = 0; trial < n_trials; ++trial) {
        auto const [u, v] = pick_2_distinct_samples_at_random(asserting_cast<tsk_id_t>(tree_sequence.num_samples()));

        auto const tskit_lca = tree_sequence.lca(u, v);

        sfkit::SampleSet samples(tree_sequence.num_samples());
        samples.add(ts_2_sf_node(0, u));
        samples.add(ts_2_sf_node(0, v));

        DAGLowestCommonAncestor lca(forest.postorder_edges());
        auto                    sfkit_lca = lca.lca(samples);

        REQUIRE(tskit_lca.size() == sfkit_lca.size());
        REQUIRE(tskit_lca.size() == forest.num_trees());
        for (TreeId tree_id = 0; tree_id < forest.num_trees(); ++tree_id) {
            CHECK(ts_2_sf_node(tree_id, asserting_cast<tsk_id_t>(tskit_lca[tree_id])) == sfkit_lca[tree_id]);
        }
    }
}
