#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>

#include "mocks/TsToSfMappingExtractor.hpp"
#include "sfkit/SuccinctForest.hpp"
#include "sfkit/assertion_levels.hpp"
#include "sfkit/dag/DAGForestCompressor.hpp"
#include "sfkit/graph/AdjacencyArrayGraph.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/graph/primitives.hpp"
#include "sfkit/utils/concepts.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace Catch::Matchers;

using sfkit::dag::DAGCompressedForest;
using sfkit::dag::DAGForestCompressor;
using sfkit::graph::NodeId;
using sfkit::samples::SampleId;
using sfkit::tskit::TSKitTree;
using sfkit::tskit::TSKitTreeSequence;

TEST_CASE("CompressedForest::is_sample() Timon", "[CompressedForest]") {
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

    std::string const   ts_file = "data/test-timon.trees";
    TSKitTreeSequence   tree_sequence(ts_file);
    DAGForestCompressor forest_compressor(tree_sequence);

    std::vector<tsk_id_t> samples = {0, 6, 10, 8, 9, 11, 12, 13, 4, 7, 14, 19, 1, 2, 18, 3, 5, 15, 16, 17};
    CHECK(samples.size() == tree_sequence.num_samples());
    CHECK(tree_sequence.num_trees() == 1);

    Ts2SfMappingExtractor ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    DAGCompressedForest   forest = forest_compressor.compress(ts_2_sf_node);
    CHECK(ts_2_sf_node.finalize_called());
    CHECK(ts_2_sf_node.process_mutations_callcnt() == 1);

    constexpr TreeId tree_id = 0;
    // Note, that the root node does not get mapped (see CompressedForest::compress() for rationale).
    for (tsk_id_t node = 0; node < asserting_cast<tsk_id_t>(tree_sequence.num_nodes() - 1); ++node) {
        auto const dag_node = ts_2_sf_node(tree_id, node);
        if (std::find(samples.begin(), samples.end(), node) != samples.end()) {
            CHECK(forest.is_sample(dag_node));
        } else {
            CHECK_FALSE(forest.is_sample(dag_node));
        }
    }
}

TEST_CASE("CompressedForest::is_sample() Scar", "[CompressedForest]") {
    //                              38
    //                    ┏━━━━━━━━━━┻━━━━━━━━━━┓
    //                   37                     ┃
    //               ┏━━━━┻━━━━┓                ┃
    //              36         ┃                ┃
    //           ┏━━━┻━━━━┓    ┃                ┃
    //           ┃        ┃   35                ┃
    //           ┃        ┃  ┏━┻━━┓             ┃
    //          34        ┃  ┃    ┃             ┃
    //      ┏━━━━┻━━━━━┓  ┃  ┃    ┃             ┃
    //      ┃          ┃  ┃  ┃    ┃            33
    //      ┃          ┃  ┃  ┃    ┃        ┏━━━━┻━━━━┓
    //      ┃          ┃  ┃  ┃    ┃       32         ┃
    //      ┃          ┃  ┃  ┃    ┃    ┏━━━┻━━━┓     ┃
    //     31          ┃  ┃  ┃    ┃    ┃       ┃     ┃
    //   ┏━━┻━━┓       ┃  ┃  ┃    ┃    ┃       ┃     ┃
    //   ┃     ┃       ┃  ┃ 30    ┃    ┃       ┃     ┃
    //   ┃     ┃       ┃  ┃ ┏┻━┓  ┃    ┃       ┃     ┃
    //  29     ┃       ┃  ┃ ┃  ┃  ┃    ┃       ┃     ┃
    // ┏━┻┓    ┃       ┃  ┃ ┃  ┃  ┃    ┃       ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃   28       ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┏━━┻━┓     ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃   27     ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃  ┏━┻━┓   ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ 26   ┃   ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┏┻━┓ ┃   ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  25     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┏┻━┓   ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  24
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┏┻━┓
    // ┃  ┃   23       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
    // ┃  ┃ ┏━━┻━┓     ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
    // ┃  ┃ ┃   22     ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
    // ┃  ┃ ┃  ┏━┻━━┓  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
    // ┃  ┃ ┃  ┃    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ 21  ┃  ┃  ┃
    // ┃  ┃ ┃  ┃    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┏┻┓ ┃  ┃  ┃
    // ┃  ┃ ┃ 20    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┃ ┃ ┃  ┃  ┃
    // ┃  ┃ ┃ ┏┻━┓  ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┃ ┃ ┃  ┃  ┃
    // 0 11 4 7 13 17 12 10 1 18 16 2 6 14 9 3 8 5 15 19

    std::string const   ts_file = "data/test-scar.trees";
    TSKitTreeSequence   tree_sequence(ts_file);
    DAGForestCompressor forest_compressor(tree_sequence);

    std::vector<tsk_id_t> samples = {0, 11, 4, 7, 13, 17, 12, 10, 1, 18, 16, 2, 6, 14, 9, 3, 8, 5, 15, 19};
    CHECK(samples.size() == tree_sequence.num_samples());
    CHECK(tree_sequence.num_trees() == 1);

    Ts2SfMappingExtractor ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    DAGCompressedForest   forest = forest_compressor.compress(ts_2_sf_node);
    CHECK(ts_2_sf_node.finalize_called());
    CHECK(ts_2_sf_node.process_mutations_callcnt() == 1);

    constexpr TreeId tree_id = 0;
    // Note, that the root node does not get mapped (see CompressedForest::compress() for rationale).
    for (tsk_id_t node = 0; node < asserting_cast<tsk_id_t>(tree_sequence.num_nodes() - 1); ++node) {
        auto const dag_node = ts_2_sf_node(tree_id, node);
        if (std::find(samples.begin(), samples.end(), node) != samples.end()) {
            CHECK(forest.is_sample(dag_node));
        } else {
            CHECK_FALSE(forest.is_sample(dag_node));
        }
    }
}

TEST_CASE("CompressedForest::is_sample() Ed", "[CompressedForest]") {
    // Too large to print here

    std::string const   ts_file = "data/test-ed.trees";
    TSKitTreeSequence   tree_sequence(ts_file);
    DAGForestCompressor forest_compressor(tree_sequence);

    std::vector<tsk_id_t> samples = {
        0,   8,   18,  148, 28,  59,  112, 78,  119, 199, 179, 190, 151, 159, 5,   32,  79,  84,  109, 182,
        154, 43,  81,  149, 195, 83,  129, 156, 6,   187, 88,  157, 181, 152, 12,  35,  114, 175, 197, 33,
        117, 15,  30,  176, 80,  137, 23,  142, 130, 172, 191, 21,  47,  107, 161, 113, 31,  174, 183, 13,
        103, 85,  122, 19,  24,  110, 69,  25,  95,  127, 89,  72,  105, 10,  111, 146, 101, 14,  11,  27,
        184, 60,  97,  141, 65,  139, 3,   50,  188, 16,  53,  133, 76,  67,  41,  51,  128, 94,  169, 77,
        82,  162, 123, 4,   44,  93,  75,  178, 9,   166, 106, 131, 40,  115, 61,  73,  180, 90,  150, 163,
        132, 124, 170, 158, 22,  57,  165, 71,  92,  34,  143, 147, 177, 192, 62,  168, 1,   104, 116, 7,
        145, 189, 42,  98,  86,  136, 155, 193, 17,  64,  37,  54,  63,  153, 118, 120, 58,  196, 48,  99,
        20,  26,  49,  121, 39,  29,  45,  144, 108, 2,   171, 38,  56,  66,  173, 100, 36,  126, 74,  46,
        55,  52,  102, 135, 185, 138, 167, 186, 70,  125, 164, 91,  198, 140, 194, 134, 160, 87,  68,  96};
    CHECK(samples.size() == tree_sequence.num_samples());
    CHECK(tree_sequence.num_trees() == 2);

    Ts2SfMappingExtractor ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    DAGCompressedForest   forest = forest_compressor.compress(ts_2_sf_node);
    CHECK(ts_2_sf_node.finalize_called());
    CHECK(ts_2_sf_node.process_mutations_callcnt() == 2);

    TSKitTree ts_tree(tree_sequence);
    TreeId    tree_id = 0;
    for (ts_tree.first(); ts_tree.is_tree(); ts_tree.next()) {
        for (auto const node: ts_tree.postorder()) {
            // Note, that the root node does not get mapped (see CompressedForest::compress() for rationale).
            if (!ts_tree.is_root(node)) {
                auto const dag_node = ts_2_sf_node(tree_id, node);
                if (std::find(samples.begin(), samples.end(), node) != samples.end()) {
                    CHECK(forest.is_sample(dag_node));
                } else {
                    CHECK_FALSE(forest.is_sample(dag_node));
                }
            }
        }
        tree_id++;
    }
}

TEST_CASE("CompressedForest TS to SF Node Mapping Timon", "[CompressedForest]") {
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

    std::string const   ts_file = "data/test-timon.trees";
    TSKitTreeSequence   tree_sequence(ts_file);
    DAGForestCompressor forest_compressor(tree_sequence);

    std::vector<tsk_id_t> samples = {0, 6, 10, 8, 9, 11, 12, 13, 4, 7, 14, 19, 1, 2, 18, 3, 5, 15, 16, 17};
    CHECK(samples.size() == tree_sequence.num_samples());
    CHECK(tree_sequence.num_trees() == 1);

    Ts2SfMappingExtractor ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    DAGCompressedForest   forest = forest_compressor.compress(ts_2_sf_node);
    CHECK(ts_2_sf_node.finalize_called());
    CHECK(ts_2_sf_node.process_mutations_callcnt() == 1);

    std::vector<NodeId> dag_samples = ts_2_sf_node(0, samples);
    CHECK(dag_samples.size() == samples.size());

    CHECK(forest.num_nodes() == 39);
    CHECK(forest.num_trees() == 1);
    CHECK(forest.num_samples() == samples.size());

    CHECK(forest.postorder_edges().num_nodes() == forest.num_nodes());
    CHECK(forest.postorder_edges().num_trees() == forest.num_trees());

    // The mapping from tree sequence nodes to forest nodes is unique if there's only one tree.
    // Note, that the root node does not get mapped (see CompressedForest::compress() for rationale).
    std::vector<tsk_id_t> all_nodes(forest.num_nodes() - 1);
    std::iota(all_nodes.begin(), all_nodes.end(), 0);
    std::vector<NodeId> dag_all_nodes = ts_2_sf_node(0, all_nodes);
    CHECK(dag_all_nodes.size() == all_nodes.size());

    auto const num_nodes = dag_all_nodes.size();
    std::sort(dag_all_nodes.begin(), dag_all_nodes.end());
    dag_all_nodes.erase(std::unique(dag_all_nodes.begin(), dag_all_nodes.end()), dag_all_nodes.end());
    CHECK(dag_all_nodes.size() == num_nodes);

    // Note, that the root node does not get mapped (see CompressedForest::compress() for rationale).
    tsk_id_t num_nodes_wo_root = asserting_cast<tsk_id_t>(tree_sequence.num_nodes() - 1);
    for (tsk_id_t ts_node = 0; ts_node < num_nodes_wo_root; ++ts_node) {
        auto const dag_node = ts_2_sf_node(0, ts_node);
        CHECK(dag_node < tree_sequence.num_nodes() - 1);
    }
}

TEST_CASE("CompressedForest TS to SF Node Mapping Shenzi", "[CompressedForest]") {
    //                              38
    //                    ┏━━━━━━━━━━┻━━━━━━━━━━┓
    //                   37                     ┃
    //               ┏━━━━┻━━━━┓                ┃
    //              36         ┃                ┃
    //           ┏━━━┻━━━━┓    ┃                ┃
    //           ┃        ┃   35                ┃
    //           ┃        ┃  ┏━┻━━┓             ┃
    //          34        ┃  ┃    ┃             ┃
    //      ┏━━━━┻━━━━━┓  ┃  ┃    ┃             ┃
    //      ┃          ┃  ┃  ┃    ┃            33
    //      ┃          ┃  ┃  ┃    ┃        ┏━━━━┻━━━━┓
    //      ┃          ┃  ┃  ┃    ┃       32         ┃
    //      ┃          ┃  ┃  ┃    ┃    ┏━━━┻━━━┓     ┃
    //     31          ┃  ┃  ┃    ┃    ┃       ┃     ┃
    //   ┏━━┻━━┓       ┃  ┃  ┃    ┃    ┃       ┃     ┃
    //   ┃     ┃       ┃  ┃ 30    ┃    ┃       ┃     ┃
    //   ┃     ┃       ┃  ┃ ┏┻━┓  ┃    ┃       ┃     ┃
    //  29     ┃       ┃  ┃ ┃  ┃  ┃    ┃       ┃     ┃
    // ┏━┻┓    ┃       ┃  ┃ ┃  ┃  ┃    ┃       ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃   28       ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┏━━┻━┓     ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃   27     ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃  ┏━┻━┓   ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ 26   ┃   ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┏┻━┓ ┃   ┃     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  25     ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┏┻━┓   ┃
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  24
    // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┏┻━┓
    // ┃  ┃   23       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
    // ┃  ┃ ┏━━┻━┓     ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
    // ┃  ┃ ┃   22     ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
    // ┃  ┃ ┃  ┏━┻━━┓  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
    // ┃  ┃ ┃  ┃    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ 21  ┃  ┃  ┃
    // ┃  ┃ ┃  ┃    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┏┻┓ ┃  ┃  ┃
    // ┃  ┃ ┃ 20    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┃ ┃ ┃  ┃  ┃
    // ┃  ┃ ┃ ┏┻━┓  ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┃ ┃ ┃  ┃  ┃
    // 0 11 4 7 13 17 12 10 1 18 16 2 6 14 9 3 8 5 15 19

    std::string const   ts_file = "data/test-shenzi.trees";
    TSKitTreeSequence   tree_sequence(ts_file);
    DAGForestCompressor forest_compressor(tree_sequence);

    Ts2SfMappingExtractor ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    DAGCompressedForest   forest = forest_compressor.compress(ts_2_sf_node);
    CHECK(ts_2_sf_node.finalize_called());
    CHECK(ts_2_sf_node.process_mutations_callcnt() == 1);

    std::vector<tsk_id_t> samples     = {0, 11, 4, 7, 13, 17, 12, 10, 1, 18, 16, 2, 6, 14, 9, 3, 8, 5, 15, 19};
    std::vector<NodeId>   dag_samples = ts_2_sf_node(0, samples);
    CHECK(dag_samples.size() == samples.size());

    CHECK(forest.num_nodes() == 39);
    CHECK(forest.num_trees() == 1);
    CHECK(forest.num_samples() == samples.size());

    CHECK(forest.postorder_edges().num_nodes() == forest.num_nodes());
    CHECK(forest.postorder_edges().num_trees() == forest.num_trees());

    // The mapping from tree sequence nodes to forest nodes is unique if there's only one tree.
    // Note, that the root node does not get mapped (see CompressedForest::compress() for rationale).
    std::vector<tsk_id_t> all_nodes(forest.num_nodes() - 1);
    std::iota(all_nodes.begin(), all_nodes.end(), 0);
    std::vector<NodeId> dag_all_nodes = ts_2_sf_node(0, all_nodes);
    CHECK(dag_all_nodes.size() == all_nodes.size());

    auto const num_nodes = dag_all_nodes.size();
    std::sort(dag_all_nodes.begin(), dag_all_nodes.end());
    dag_all_nodes.erase(std::unique(dag_all_nodes.begin(), dag_all_nodes.end()), dag_all_nodes.end());
    CHECK(dag_all_nodes.size() == num_nodes);

    // Note, that the root node does not get mapped (see CompressedForest::compress() for rationale).
    tsk_id_t num_nodes_wo_root = asserting_cast<tsk_id_t>(tree_sequence.num_nodes() - 1);
    for (tsk_id_t ts_node = 0; ts_node < num_nodes_wo_root; ++ts_node) {
        auto const dag_node = ts_2_sf_node(0, ts_node);
        CHECK(dag_node < tree_sequence.num_nodes() - 1);
    }
}

TEST_CASE("Subtrees are created only once Example I", "[CompresedForest]") {
    /*
    0.25┊     8   ┊         ┊         ┊
        ┊   ┏━┻━┓ ┊         ┊         ┊
    0.20┊   ┃   ┃ ┊         ┊   7     ┊
        ┊   ┃   ┃ ┊         ┊ ┏━┻━┓   ┊
    0.17┊   6   ┃ ┊   6     ┊ ┃   ┃   ┊
        ┊ ┏━┻┓  ┃ ┊ ┏━┻━┓   ┊ ┃   ┃   ┊
    0.09┊ ┃  5  ┃ ┊ ┃   5   ┊ ┃   5   ┊
        ┊ ┃ ┏┻┓ ┃ ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊
    0.07┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊
        ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊
    0.00┊ 0 1 3 2 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊
    */

    tsk_treeseq_t tree_sequence;

    tsk_treeseq_from_text(
        &tree_sequence,
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

    sfkit::DAGSuccinctForestNumeric forest(tree_sequence); // Takes ownership

    CHECK(forest.num_trees() == 3);
    CHECK(forest.num_samples() == 4);
    // Even though the last two trees are identical, the root exists twice (see forest compression).
    CHECK(forest.num_unique_subtrees() == 11);
}
