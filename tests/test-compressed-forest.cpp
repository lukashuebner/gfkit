#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>

#include "tdt/assertion_levels.hpp"
#include "tdt/graph/adjacency-array-graph.hpp"
#include "tdt/graph/common.hpp"
#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/load/forest-compressor.hpp"
#include "tdt/utils/concepts.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace Catch::Matchers;

// TEST_CASE("CompressedForest::is_sample() I", "[CompressedForest]") {
//     //                                      38
//     //                                ┏━━━━━━┻━━━━━━━┓
//     //                               37              ┃
//     //                         ┏━━━━━━┻━━━━━━┓       ┃
//     //                        36             ┃       ┃
//     //                 ┏━━━━━━━┻━━━━━━━━┓    ┃       ┃
//     //                 ┃                ┃   35       ┃
//     //                 ┃                ┃   ┏┻┓      ┃
//     //                 ┃                ┃   ┃ ┃     34
//     //                 ┃                ┃   ┃ ┃    ┏━┻━┓
//     //                33                ┃   ┃ ┃    ┃   ┃
//     //        ┏━━━━━━━━┻━━━━━━━━┓       ┃   ┃ ┃    ┃   ┃
//     //       32                 ┃       ┃   ┃ ┃    ┃   ┃
//     //    ┏━━━┻━━━┓             ┃       ┃   ┃ ┃    ┃   ┃
//     //   31       ┃             ┃       ┃   ┃ ┃    ┃   ┃
//     //  ┏━┻━┓     ┃             ┃       ┃   ┃ ┃    ┃   ┃
//     //  ┃   ┃     ┃             ┃       ┃   ┃ ┃   30   ┃
//     //  ┃   ┃     ┃             ┃       ┃   ┃ ┃  ┏━┻┓  ┃
//     //  ┃   ┃    29             ┃       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┏━━━┻━━━┓         ┃       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃      28         ┃       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┏━━┻━━┓      ┃       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┃     ┃     27       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┃     ┃   ┏━━┻━━┓    ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┃     ┃  26     ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┃     ┃ ┏━┻━┓   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃   25     ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃  ┏━┻━━┓  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // 24   ┃ ┃  ┃    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // ┏┻┓  ┃ ┃  ┃    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ 23    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┏┻━┓  ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃   ┃   ┃   22   ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃   ┃   ┃  ┏━┻━┓ ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃  21   ┃  ┃   ┃ ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┏━┻┓  ┃  ┃   ┃ ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┃  ┃  ┃ 20   ┃ ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┃  ┃  ┃ ┏┻┓  ┃ ┃ ┃  ┃  ┃  ┃
//     // 0 6 10 8 9 11 12 13 4 7 14 19 1 2 18 3 5 15 16 17

//     std::string const ts_file = "data/compressed-forest-is-sample.trees";
//     TSKitTreeSequence tree_sequence(ts_file);
//     ForestCompressor  forest_compressor(tree_sequence);
//     CompressedForest  forest = forest_compressor.compress();

//     std::vector<tsk_id_t> samples = {0, 6, 10, 8, 9, 11, 12, 13, 4, 7, 14, 19, 1, 2, 18, 3, 5, 15, 16, 17};
//     CHECK(samples.size() == tree_sequence.num_samples());
//     CHECK(tree_sequence.num_trees() == 1);
//     for (tsk_id_t node = 0; node <= 38; ++node) {
//         SubtreeId dag_node = forest_compressor.ts_node2cf_subtree(0, node);
//         if (std::find(samples.begin(), samples.end(), node) != samples.end()) {
//             CHECK(forest.is_sample(dag_node));
//         } else {
//             CHECK_FALSE(forest.is_sample(dag_node));
//         }
//     }
// }

// TEST_CASE("CompressedForest::is_sample() II", "[CompressedForest]") {
//     //                              38
//     //                    ┏━━━━━━━━━━┻━━━━━━━━━━┓
//     //                   37                     ┃
//     //               ┏━━━━┻━━━━┓                ┃
//     //              36         ┃                ┃
//     //           ┏━━━┻━━━━┓    ┃                ┃
//     //           ┃        ┃   35                ┃
//     //           ┃        ┃  ┏━┻━━┓             ┃
//     //          34        ┃  ┃    ┃             ┃
//     //      ┏━━━━┻━━━━━┓  ┃  ┃    ┃             ┃
//     //      ┃          ┃  ┃  ┃    ┃            33
//     //      ┃          ┃  ┃  ┃    ┃        ┏━━━━┻━━━━┓
//     //      ┃          ┃  ┃  ┃    ┃       32         ┃
//     //      ┃          ┃  ┃  ┃    ┃    ┏━━━┻━━━┓     ┃
//     //     31          ┃  ┃  ┃    ┃    ┃       ┃     ┃
//     //   ┏━━┻━━┓       ┃  ┃  ┃    ┃    ┃       ┃     ┃
//     //   ┃     ┃       ┃  ┃ 30    ┃    ┃       ┃     ┃
//     //   ┃     ┃       ┃  ┃ ┏┻━┓  ┃    ┃       ┃     ┃
//     //  29     ┃       ┃  ┃ ┃  ┃  ┃    ┃       ┃     ┃
//     // ┏━┻┓    ┃       ┃  ┃ ┃  ┃  ┃    ┃       ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃   28       ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┏━━┻━┓     ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃   27     ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃  ┏━┻━┓   ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ 26   ┃   ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┏┻━┓ ┃   ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  25     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┏┻━┓   ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  24
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┏┻━┓
//     // ┃  ┃   23       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
//     // ┃  ┃ ┏━━┻━┓     ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
//     // ┃  ┃ ┃   22     ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
//     // ┃  ┃ ┃  ┏━┻━━┓  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
//     // ┃  ┃ ┃  ┃    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ 21  ┃  ┃  ┃
//     // ┃  ┃ ┃  ┃    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┏┻┓ ┃  ┃  ┃
//     // ┃  ┃ ┃ 20    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┃ ┃ ┃  ┃  ┃
//     // ┃  ┃ ┃ ┏┻━┓  ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┃ ┃ ┃  ┃  ┃
//     // 0 11 4 7 13 17 12 10 1 18 16 2 6 14 9 3 8 5 15 19

//     std::string const ts_file = "data/allele-frequency-spectrum-simple-example-2.trees";
//     TSKitTreeSequence tree_sequence(ts_file);
//     ForestCompressor  forest_compressor(tree_sequence);
//     CompressedForest  forest = forest_compressor.compress();

//     std::vector<tsk_id_t> samples = {0, 11, 4, 7, 13, 17, 12, 10, 1, 18, 16, 2, 6, 14, 9, 3, 8, 5, 15, 19};
//     CHECK(samples.size() == tree_sequence.num_samples());
//     CHECK(tree_sequence.num_trees() == 1);
//     for (tsk_id_t node = 0; node <= 38; ++node) {
//         SubtreeId dag_node = forest_compressor.ts_node2cf_subtree(0, node);
//         if (std::find(samples.begin(), samples.end(), node) != samples.end()) {
//             CHECK(forest.is_sample(dag_node));
//         } else {
//             CHECK_FALSE(forest.is_sample(dag_node));
//         }
//     }
// }

// TEST_CASE("CompressedForest::is_sample() III", "[CompressedForest]") {
//     // Too large to print here

//     std::string const     ts_file = "data/allele-frequency-spectrum-simple-example-4.trees";
//     TSKitTreeSequence     tree_sequence(ts_file);
//     ForestCompressor      forest_compressor(tree_sequence);
//     CompressedForest      forest  = forest_compressor.compress();
//     std::vector<tsk_id_t> samples = {
//         0,   8,   18,  148, 28,  59,  112, 78,  119, 199, 179, 190, 151, 159, 5,   32,  79,  84,  109, 182,
//         154, 43,  81,  149, 195, 83,  129, 156, 6,   187, 88,  157, 181, 152, 12,  35,  114, 175, 197, 33,
//         117, 15,  30,  176, 80,  137, 23,  142, 130, 172, 191, 21,  47,  107, 161, 113, 31,  174, 183, 13,
//         103, 85,  122, 19,  24,  110, 69,  25,  95,  127, 89,  72,  105, 10,  111, 146, 101, 14,  11,  27,
//         184, 60,  97,  141, 65,  139, 3,   50,  188, 16,  53,  133, 76,  67,  41,  51,  128, 94,  169, 77,
//         82,  162, 123, 4,   44,  93,  75,  178, 9,   166, 106, 131, 40,  115, 61,  73,  180, 90,  150, 163,
//         132, 124, 170, 158, 22,  57,  165, 71,  92,  34,  143, 147, 177, 192, 62,  168, 1,   104, 116, 7,
//         145, 189, 42,  98,  86,  136, 155, 193, 17,  64,  37,  54,  63,  153, 118, 120, 58,  196, 48,  99,
//         20,  26,  49,  121, 39,  29,  45,  144, 108, 2,   171, 38,  56,  66,  173, 100, 36,  126, 74,  46,
//         55,  52,  102, 135, 185, 138, 167, 186, 70,  125, 164, 91,  198, 140, 194, 134, 160, 87,  68,  96};
//     CHECK(samples.size() == tree_sequence.num_samples());
//     CHECK(tree_sequence.num_trees() == 2);
//     for (tsk_id_t node = 0; node <= 38; ++node) {
//         SubtreeId dag_node = forest_compressor.ts_node2cf_subtree(0, node);
//         if (std::find(samples.begin(), samples.end(), node) != samples.end()) {
//             CHECK(forest.is_sample(dag_node));
//         } else {
//             CHECK_FALSE(forest.is_sample(dag_node));
//         }
//     }
// }

// TEST_CASE("CompressedForest::{num_,}samples_below() I", "[CompressedForest]") {
//     //                                      38
//     //                                ┏━━━━━━┻━━━━━━━┓
//     //                               37              ┃
//     //                         ┏━━━━━━┻━━━━━━┓       ┃
//     //                        36             ┃       ┃
//     //                 ┏━━━━━━━┻━━━━━━━━┓    ┃       ┃
//     //                 ┃                ┃   35       ┃
//     //                 ┃                ┃   ┏┻┓      ┃
//     //                 ┃                ┃   ┃ ┃     34
//     //                 ┃                ┃   ┃ ┃    ┏━┻━┓
//     //                33                ┃   ┃ ┃    ┃   ┃
//     //        ┏━━━━━━━━┻━━━━━━━━┓       ┃   ┃ ┃    ┃   ┃
//     //       32                 ┃       ┃   ┃ ┃    ┃   ┃
//     //    ┏━━━┻━━━┓             ┃       ┃   ┃ ┃    ┃   ┃
//     //   31       ┃             ┃       ┃   ┃ ┃    ┃   ┃
//     //  ┏━┻━┓     ┃             ┃       ┃   ┃ ┃    ┃   ┃
//     //  ┃   ┃     ┃             ┃       ┃   ┃ ┃   30   ┃
//     //  ┃   ┃     ┃             ┃       ┃   ┃ ┃  ┏━┻┓  ┃
//     //  ┃   ┃    29             ┃       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┏━━━┻━━━┓         ┃       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃      28         ┃       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┏━━┻━━┓      ┃       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┃     ┃     27       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┃     ┃   ┏━━┻━━┓    ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┃     ┃  26     ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┃     ┃ ┏━┻━┓   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃   25     ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃  ┏━┻━━┓  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // 24   ┃ ┃  ┃    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // ┏┻┓  ┃ ┃  ┃    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ 23    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┏┻━┓  ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃   ┃   ┃   22   ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃   ┃   ┃  ┏━┻━┓ ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃  21   ┃  ┃   ┃ ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┏━┻┓  ┃  ┃   ┃ ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┃  ┃  ┃ 20   ┃ ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┃  ┃  ┃ ┏┻┓  ┃ ┃ ┃  ┃  ┃  ┃
//     // 0 6 10 8 9 11 12 13 4 7 14 19 1 2 18 3 5 15 16 17

//     std::string const ts_file = "data/compressed-forest-is-sample.trees";
//     TSKitTreeSequence tree_sequence(ts_file);
//     ForestCompressor  forest_compressor(tree_sequence);
//     CompressedForest  forest = forest_compressor.compress();

//     forest.compute_num_samples_below();
//     // forest.compute_samples_below();

//     std::vector<tsk_id_t> samples     = {0, 6, 10, 8, 9, 11, 12, 13, 4, 7, 14, 19, 1, 2, 18, 3, 5, 15, 16, 17};
//     std::vector<NodeId>   dag_samples = forest_compressor.ts_node2cf_subtree(0, samples);

//     for (NodeId dag_sample: dag_samples) {
//         CHECK(forest.is_sample(dag_sample));
//         // CHECK_THAT(forest.samples_below(dag_sample), UnorderedRangeEquals(std::vector<NodeId>{dag_sample}));
//     }

//     CHECK(samples.size() == tree_sequence.num_samples());
//     CHECK(samples.size() == dag_samples.size());
//     CHECK(tree_sequence.num_trees() == 1);
//     for (tsk_id_t node = 0; node <= 38; ++node) {
//         SubtreeId dag_node = forest_compressor.ts_node2cf_subtree(0, node);
//         if (std::find(samples.begin(), samples.end(), node) != samples.end()) {
//             CHECK(forest.is_sample(dag_node));
//             // CHECK(forest.samples_below(dag_node).size() == 1);
//         } else {
//             CHECK_FALSE(forest.is_sample(dag_node));
//             // CHECK(forest.samples_below(dag_node).size() > 1);
//         }
//     }

//     // CHECK_THAT(forest.samples_below(forest_compressor.ts_node2cf_subtree(0, 38)), UnorderedRangeEquals(dag_samples));
//     // CHECK_THAT(
//     //     forest.samples_below(forest_compressor.ts_node2cf_subtree(0, 20)),
//     //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(0, std::vector<NodeId>{1, 2}))
//     // );
//     // CHECK_THAT(
//     //     forest.samples_below(forest_compressor.ts_node2cf_subtree(0, 27)),
//     //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(0, std::vector<NodeId>{4, 7, 14, 19}))
//     // );
//     // CHECK_THAT(
//     //     forest.samples_below(forest_compressor.ts_node2cf_subtree(0, 29)),
//     //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(0, std::vector<NodeId>{8, 9, 11, 12, 13}))
//     // );
//     //
//     // for (tsk_id_t ts_node = 0; ts_node <= 38; ++ts_node) {
//     //     auto const dag_node = forest_compressor.ts_node2cf_subtree(0, ts_node);
//     //     CHECK(forest.samples_below(dag_node).size() == forest.num_samples_below(dag_node));
//     // }
// }

// TEST_CASE("CompressedForest::{num_,}samples_below() II", "[CompressedForest]") {
//     //                              38
//     //                    ┏━━━━━━━━━━┻━━━━━━━━━━┓
//     //                   37                     ┃
//     //               ┏━━━━┻━━━━┓                ┃
//     //              36         ┃                ┃
//     //           ┏━━━┻━━━━┓    ┃                ┃
//     //           ┃        ┃   35                ┃
//     //           ┃        ┃  ┏━┻━━┓             ┃
//     //          34        ┃  ┃    ┃             ┃
//     //      ┏━━━━┻━━━━━┓  ┃  ┃    ┃             ┃
//     //      ┃          ┃  ┃  ┃    ┃            33
//     //      ┃          ┃  ┃  ┃    ┃        ┏━━━━┻━━━━┓
//     //      ┃          ┃  ┃  ┃    ┃       32         ┃
//     //      ┃          ┃  ┃  ┃    ┃    ┏━━━┻━━━┓     ┃
//     //     31          ┃  ┃  ┃    ┃    ┃       ┃     ┃
//     //   ┏━━┻━━┓       ┃  ┃  ┃    ┃    ┃       ┃     ┃
//     //   ┃     ┃       ┃  ┃ 30    ┃    ┃       ┃     ┃
//     //   ┃     ┃       ┃  ┃ ┏┻━┓  ┃    ┃       ┃     ┃
//     //  29     ┃       ┃  ┃ ┃  ┃  ┃    ┃       ┃     ┃
//     // ┏━┻┓    ┃       ┃  ┃ ┃  ┃  ┃    ┃       ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃   28       ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┏━━┻━┓     ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃   27     ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃  ┏━┻━┓   ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ 26   ┃   ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┏┻━┓ ┃   ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  25     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┏┻━┓   ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  24
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┏┻━┓
//     // ┃  ┃   23       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
//     // ┃  ┃ ┏━━┻━┓     ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
//     // ┃  ┃ ┃   22     ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
//     // ┃  ┃ ┃  ┏━┻━━┓  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
//     // ┃  ┃ ┃  ┃    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ 21  ┃  ┃  ┃
//     // ┃  ┃ ┃  ┃    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┏┻┓ ┃  ┃  ┃
//     // ┃  ┃ ┃ 20    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┃ ┃ ┃  ┃  ┃
//     // ┃  ┃ ┃ ┏┻━┓  ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┃ ┃ ┃  ┃  ┃
//     // 0 11 4 7 13 17 12 10 1 18 16 2 6 14 9 3 8 5 15 19

//     std::string const ts_file = "data/allele-frequency-spectrum-simple-example-2.trees";
//     TSKitTreeSequence tree_sequence(ts_file);
//     ForestCompressor  forest_compressor(tree_sequence);
//     CompressedForest  forest = forest_compressor.compress();
//     forest.compute_num_samples_below();
//     // forest.compute_samples_below();

//     std::vector<tsk_id_t> samples     = {0, 11, 4, 7, 13, 17, 12, 10, 1, 18, 16, 2, 6, 14, 9, 3, 8, 5, 15, 19};
//     std::vector<NodeId>   dag_samples = forest_compressor.ts_node2cf_subtree(0, samples);

//     // CHECK_THAT(forest.samples_below(38), UnorderedRangeEquals(dag_samples));
//     for (NodeId dag_sample: dag_samples) {
//         CHECK(forest.is_sample(dag_sample));
//         // CHECK_THAT(forest.samples_below(dag_sample), UnorderedRangeEquals(std::vector<NodeId>{dag_sample}));
//     }
//     // CHECK_THAT(
//     //     forest.samples_below(forest_compressor.ts_node2cf_subtree(0, 33)),
//     //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(0, std::vector<NodeId>{2, 6, 14, 9, 3, 8, 5, 15, 19}))
//     // );

//     // for (tsk_id_t ts_node = 0; ts_node <= 38; ++ts_node) {
//     //     auto const dag_node = forest_compressor.ts_node2cf_subtree(0, ts_node);
//     //     CHECK(forest.samples_below(dag_node).size() == forest.num_samples_below(dag_node));
//     // }
// }

// TEST_CASE("CompressedForest::{num_,}samples_below() III", "[CompressedForest]") {
//     //                                    39
//     //                         ┏━━━━━━━━━━━┻━━━━━━━━━━━┓
//     //                        38                       ┃
//     //           ┏━━━━━━━━━━━━━┻━━━━━━━━━━━━━━┓        ┃
//     //           ┃                           37        ┃
//     //           ┃                       ┏━━━━┻━━━━┓   ┃
//     //          36                       ┃         ┃   ┃
//     //      ┏━━━━┻━━━┓                   ┃         ┃   ┃
//     //      ┃       35                   ┃         ┃   ┃
//     //      ┃    ┏━━━┻━━━┓               ┃         ┃   ┃
//     //     34    ┃       ┃               ┃         ┃   ┃
//     //    ┏━┻━┓  ┃       ┃               ┃         ┃   ┃
//     //    ┃   ┃  ┃       ┃              33         ┃   ┃
//     //    ┃   ┃  ┃       ┃          ┏━━━━┻━━━━━┓   ┃   ┃
//     //    ┃   ┃  ┃      32          ┃          ┃   ┃   ┃
//     //    ┃   ┃  ┃    ┏━━┻━┓        ┃          ┃   ┃   ┃
//     //    ┃   ┃  ┃    ┃   31        ┃          ┃   ┃   ┃
//     //    ┃   ┃  ┃    ┃  ┏━┻┓       ┃          ┃   ┃   ┃
//     //    ┃   ┃  ┃    ┃  ┃  ┃       ┃          ┃  30   ┃
//     //    ┃   ┃  ┃    ┃  ┃  ┃       ┃          ┃  ┏┻━┓ ┃
//     //    ┃   ┃  ┃    ┃  ┃  ┃       ┃         29  ┃  ┃ ┃
//     //    ┃   ┃  ┃    ┃  ┃  ┃       ┃        ┏━┻┓ ┃  ┃ ┃
//     //    ┃   ┃ 28    ┃  ┃  ┃       ┃        ┃  ┃ ┃  ┃ ┃
//     //    ┃   ┃ ┏┻━┓  ┃  ┃  ┃       ┃        ┃  ┃ ┃  ┃ ┃
//     //    ┃   ┃ ┃  ┃  ┃  ┃  ┃      27        ┃  ┃ ┃  ┃ ┃
//     //    ┃   ┃ ┃  ┃  ┃  ┃  ┃    ┏━━┻━━┓     ┃  ┃ ┃  ┃ ┃
//     //   26   ┃ ┃  ┃  ┃  ┃  ┃    ┃     ┃     ┃  ┃ ┃  ┃ ┃
//     //  ┏━┻━┓ ┃ ┃  ┃  ┃  ┃  ┃    ┃     ┃     ┃  ┃ ┃  ┃ ┃
//     //  ┃   ┃ ┃ ┃  ┃  ┃  ┃  ┃    ┃    25     ┃  ┃ ┃  ┃ ┃
//     //  ┃   ┃ ┃ ┃  ┃  ┃  ┃  ┃    ┃   ┏━┻━┓   ┃  ┃ ┃  ┃ ┃
//     // 23   ┃ ┃ ┃  ┃  ┃  ┃  ┃    ┃   ┃   ┃   ┃  ┃ ┃  ┃ ┃
//     // ┏┻┓  ┃ ┃ ┃  ┃  ┃  ┃  ┃    ┃   ┃   ┃   ┃  ┃ ┃  ┃ ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃  ┃   22   ┃   ┃   ┃  ┃ ┃  ┃ ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃  ┃  ┏━┻━┓ ┃   ┃   ┃  ┃ ┃  ┃ ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃  ┃  ┃   ┃ ┃  21   ┃  ┃ ┃  ┃ ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃  ┃  ┃   ┃ ┃ ┏━┻┓  ┃  ┃ ┃  ┃ ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃  ┃ 20   ┃ ┃ ┃  ┃  ┃  ┃ ┃  ┃ ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃  ┃ ┏┻┓  ┃ ┃ ┃  ┃  ┃  ┃ ┃  ┃ ┃
//     // 0 8 11 6 5 19 12 13 15 1 2 18 4 7 14 16 17 9 10 3

//     std::string const ts_file = "data/allele-frequency-spectrum-simple-example-6.trees";
//     TSKitTreeSequence tree_sequence(ts_file);
//     ForestCompressor  forest_compressor(tree_sequence);
//     CompressedForest  forest = forest_compressor.compress();

//     forest.compute_num_samples_below();
//     // forest.compute_samples_below();

//     std::vector<tsk_id_t> samples     = {0, 8, 11, 6, 5, 19, 12, 13, 15, 1, 2, 18, 4, 7, 14, 16, 17, 9, 10, 3};
//     std::vector<NodeId>   dag_samples = forest_compressor.ts_node2cf_subtree(0, samples);

//     // The samples of both trees should be mapped to the same nodes in the compressed forest DAG.
//     CHECK_THAT(dag_samples, UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(1, samples)));

//     // Each tree has its own root in the compressed forest, even if they would be exactly identical (they aren't in this
//     // example). This simplifies writing code using the compressed forest as the mapping tree <-> root in DAG is
//     // bijective.
//     tsk_id_t const ts_root = 39;
//     CHECK(forest_compressor.ts_node2cf_subtree(0, ts_root) != forest_compressor.ts_node2cf_subtree(1, ts_root));

//     SECTION("Check Tree 0") {
//         // size_t const tree0 = 0;
//         // TODO rename dag_ to cf_ everywhere
//         // NodeId const cf_root = forest_compressor.ts_node2cf_subtree(tree0, ts_root);
//         // CHECK_THAT(forest.samples_below(cf_root), UnorderedRangeEquals(dag_samples));
//         for (NodeId dag_sample: dag_samples) {
//             CHECK(forest.is_sample(dag_sample));
//             // CHECK_THAT(forest.samples_below(dag_sample), UnorderedRangeEquals(std::vector<NodeId>{dag_sample}));
//         }
//         // CHECK_THAT(
//         //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree0, 33)),
//         //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(0, std::vector<NodeId>{1, 2, 18, 4, 7, 14, 16, 17}))
//         // );
//         // CHECK_THAT(
//         //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree0, 26)),
//         //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree0, std::vector<NodeId>{0, 8, 11}))
//         // );
//         // CHECK_THAT(
//         //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree0, 32)),
//         //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree0, std::vector<NodeId>{12, 13, 15}))
//         // );
//         // CHECK_THAT(
//         //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree0, 36)),
//         //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree0, std::vector<NodeId>{0, 8, 11, 6, 5, 19, 12, 13,
//         //     15}))
//         // );

//         // for (tsk_id_t ts_node = 0; ts_node <= 39; ++ts_node) {
//         //     auto const dag_node = forest_compressor.ts_node2cf_subtree(0, ts_node);
//         //     CHECK(forest.samples_below(dag_node).size() == forest.num_samples_below(dag_node));
//         // }
//     }

//     //                                     39
//     //                           ┏━━━━━━━━━━┻━━━━━━━━━━┓
//     //                          38                     ┃
//     //           ┏━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━┓      ┃
//     //           ┃                             37      ┃
//     //           ┃                         ┏━━━━┻━━━━┓ ┃
//     //          36                         ┃         ┃ ┃
//     //      ┏━━━━┻━━━┓                     ┃         ┃ ┃
//     //      ┃       35                     ┃         ┃ ┃
//     //      ┃    ┏━━━┻━━━━┓                ┃         ┃ ┃
//     //     34    ┃        ┃                ┃         ┃ ┃
//     //    ┏━┻━┓  ┃        ┃                ┃         ┃ ┃
//     //    ┃   ┃  ┃        ┃               33         ┃ ┃
//     //    ┃   ┃  ┃        ┃           ┏━━━━┻━━━━━┓   ┃ ┃
//     //    ┃   ┃  ┃       32           ┃          ┃   ┃ ┃
//     //    ┃   ┃  ┃     ┏━━┻━━┓        ┃          ┃   ┃ ┃
//     //    ┃   ┃  ┃     ┃    31        ┃          ┃   ┃ ┃
//     //    ┃   ┃  ┃     ┃   ┏━┻┓       ┃          ┃   ┃ ┃
//     //    ┃   ┃  ┃     ┃   ┃  ┃       ┃         29   ┃ ┃
//     //    ┃   ┃  ┃     ┃   ┃  ┃       ┃        ┏━┻┓  ┃ ┃
//     //    ┃   ┃ 28     ┃   ┃  ┃       ┃        ┃  ┃  ┃ ┃
//     //    ┃   ┃ ┏┻━┓   ┃   ┃  ┃       ┃        ┃  ┃  ┃ ┃
//     //    ┃   ┃ ┃  ┃   ┃   ┃  ┃      27        ┃  ┃  ┃ ┃
//     //    ┃   ┃ ┃  ┃   ┃   ┃  ┃    ┏━━┻━━┓     ┃  ┃  ┃ ┃
//     //   26   ┃ ┃  ┃   ┃   ┃  ┃    ┃     ┃     ┃  ┃  ┃ ┃
//     //  ┏━┻━┓ ┃ ┃  ┃   ┃   ┃  ┃    ┃     ┃     ┃  ┃  ┃ ┃
//     //  ┃   ┃ ┃ ┃  ┃   ┃   ┃  ┃    ┃    25     ┃  ┃  ┃ ┃
//     //  ┃   ┃ ┃ ┃  ┃   ┃   ┃  ┃    ┃   ┏━┻━┓   ┃  ┃  ┃ ┃
//     //  ┃   ┃ ┃ ┃  ┃  24   ┃  ┃    ┃   ┃   ┃   ┃  ┃  ┃ ┃
//     //  ┃   ┃ ┃ ┃  ┃ ┏━┻┓  ┃  ┃    ┃   ┃   ┃   ┃  ┃  ┃ ┃
//     // 23   ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃    ┃   ┃   ┃   ┃  ┃  ┃ ┃
//     // ┏┻┓  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃    ┃   ┃   ┃   ┃  ┃  ┃ ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃   22   ┃   ┃   ┃  ┃  ┃ ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┏━┻━┓ ┃   ┃   ┃  ┃  ┃ ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃   ┃ ┃  21   ┃  ┃  ┃ ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃   ┃ ┃ ┏━┻┓  ┃  ┃  ┃ ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃ 20   ┃ ┃ ┃  ┃  ┃  ┃  ┃ ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃ ┏┻┓  ┃ ┃ ┃  ┃  ┃  ┃  ┃ ┃
//     // 0 8 11 6 5 19 9 12 13 15 1 2 18 4 7 14 16 17 10 3
//     SECTION("Check Tree 1") {
//         // size_t const tree1   = 1;
//         // NodeId const cf_root = forest_compressor.ts_node2cf_subtree(tree1, ts_root);
//         // CHECK_THAT(forest.samples_below(cf_root), UnorderedRangeEquals(dag_samples));
//         for (NodeId dag_sample: dag_samples) {
//             CHECK(forest.is_sample(dag_sample));
//             // CHECK_THAT(forest.samples_below(dag_sample), UnorderedRangeEquals(std::vector<NodeId>{dag_sample}));
//         }
//         // CHECK_THAT(
//         //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree1, 33)),
//         //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree1, std::vector<NodeId>{1, 2, 18, 4, 7, 14, 16, 17}))
//         // );
//         // CHECK_THAT(
//         //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree1, 26)),
//         //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree1, std::vector<NodeId>{0, 8, 11}))
//         // );
//         // CHECK_THAT(
//         //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree1, 32)),
//         //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree1, std::vector<NodeId>{9, 12, 13, 15}))
//         // );
//         // CHECK_THAT(
//         //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree1, 36)),
//         //     UnorderedRangeEquals(
//         //         forest_compressor.ts_node2cf_subtree(tree1, std::vector<NodeId>{0, 8, 11, 6, 5, 19, 9, 12, 13, 15})
//         //     )
//         // );

//         // for (tsk_id_t ts_node = 0; ts_node <= 39; ++ts_node) {
//         //     auto const dag_node = forest_compressor.ts_node2cf_subtree(0, ts_node);
//         //     CHECK(forest.samples_below(dag_node).size() == forest.num_samples_below(dag_node));
//         // }
//     }
// }

// TODO Use mocking to test the mapper passed to CompressedSequenceStorage
// TEST_CASE("CompressedForest::ts_node2cf_subtree() I", "[CompressedForest]") {
//     //                                      38
//     //                                ┏━━━━━━┻━━━━━━━┓
//     //                               37              ┃
//     //                         ┏━━━━━━┻━━━━━━┓       ┃
//     //                        36             ┃       ┃
//     //                 ┏━━━━━━━┻━━━━━━━━┓    ┃       ┃
//     //                 ┃                ┃   35       ┃
//     //                 ┃                ┃   ┏┻┓      ┃
//     //                 ┃                ┃   ┃ ┃     34
//     //                 ┃                ┃   ┃ ┃    ┏━┻━┓
//     //                33                ┃   ┃ ┃    ┃   ┃
//     //        ┏━━━━━━━━┻━━━━━━━━┓       ┃   ┃ ┃    ┃   ┃
//     //       32                 ┃       ┃   ┃ ┃    ┃   ┃
//     //    ┏━━━┻━━━┓             ┃       ┃   ┃ ┃    ┃   ┃
//     //   31       ┃             ┃       ┃   ┃ ┃    ┃   ┃
//     //  ┏━┻━┓     ┃             ┃       ┃   ┃ ┃    ┃   ┃
//     //  ┃   ┃     ┃             ┃       ┃   ┃ ┃   30   ┃
//     //  ┃   ┃     ┃             ┃       ┃   ┃ ┃  ┏━┻┓  ┃
//     //  ┃   ┃    29             ┃       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┏━━━┻━━━┓         ┃       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃      28         ┃       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┏━━┻━━┓      ┃       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┃     ┃     27       ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┃     ┃   ┏━━┻━━┓    ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┃     ┃  26     ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃    ┃     ┃ ┏━┻━┓   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃   25     ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     //  ┃   ┃ ┃  ┏━┻━━┓  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // 24   ┃ ┃  ┃    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // ┏┻┓  ┃ ┃  ┃    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ 23    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┏┻━┓  ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃   ┃   ┃   22   ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃   ┃   ┃  ┏━┻━┓ ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃  21   ┃  ┃   ┃ ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┏━┻┓  ┃  ┃   ┃ ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┃  ┃  ┃ 20   ┃ ┃ ┃  ┃  ┃  ┃
//     // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┃  ┃  ┃ ┏┻┓  ┃ ┃ ┃  ┃  ┃  ┃
//     // 0 6 10 8 9 11 12 13 4 7 14 19 1 2 18 3 5 15 16 17

//     std::string const ts_file = "data/compressed-forest-is-sample.trees";
//     TSKitTreeSequence tree_sequence(ts_file);
//     ForestCompressor  forest_compressor(tree_sequence);
//     CompressedForest  forest = forest_compressor.compress();

//     std::vector<tsk_id_t> samples     = {0, 6, 10, 8, 9, 11, 12, 13, 4, 7, 14, 19, 1, 2, 18, 3, 5, 15, 16, 17};
//     std::vector<NodeId>   dag_samples = forest_compressor.ts_node2cf_subtree(0, samples);
//     CHECK(dag_samples.size() == samples.size());

//     CHECK(forest.num_nodes() == 39);
//     CHECK(forest.num_trees() == 1);
//     CHECK(forest.num_samples() == samples.size());

//     CHECK(forest.postorder_edges().num_nodes() == forest.num_nodes());
//     CHECK(forest.postorder_edges().num_trees() == forest.num_trees());

//     // The mapping from tree sequence nodes to forest nodes is unique if there's only one tree.
//     std::vector<tsk_id_t> all_nodes(forest.num_nodes());
//     std::iota(all_nodes.begin(), all_nodes.end(), 0);
//     std::vector<NodeId> dag_all_nodes = forest_compressor.ts_node2cf_subtree(0, all_nodes);
//     CHECK(dag_all_nodes.size() == all_nodes.size());

//     auto const num_nodes = dag_all_nodes.size();
//     std::sort(dag_all_nodes.begin(), dag_all_nodes.end());
//     dag_all_nodes.erase(std::unique(dag_all_nodes.begin(), dag_all_nodes.end()), dag_all_nodes.end());
//     CHECK(dag_all_nodes.size() == num_nodes);

//     for (tsk_id_t ts_node = 0; ts_node <= 38; ++ts_node) {
//         auto const dag_node = forest_compressor.ts_node2cf_subtree(0, ts_node);
//         CHECK(dag_node <= 38);
//     }
// }

// TEST_CASE("CompressedForest::ts_node2cf_subtree() II", "[CompressedForest]") {
//     //                              38
//     //                    ┏━━━━━━━━━━┻━━━━━━━━━━┓
//     //                   37                     ┃
//     //               ┏━━━━┻━━━━┓                ┃
//     //              36         ┃                ┃
//     //           ┏━━━┻━━━━┓    ┃                ┃
//     //           ┃        ┃   35                ┃
//     //           ┃        ┃  ┏━┻━━┓             ┃
//     //          34        ┃  ┃    ┃             ┃
//     //      ┏━━━━┻━━━━━┓  ┃  ┃    ┃             ┃
//     //      ┃          ┃  ┃  ┃    ┃            33
//     //      ┃          ┃  ┃  ┃    ┃        ┏━━━━┻━━━━┓
//     //      ┃          ┃  ┃  ┃    ┃       32         ┃
//     //      ┃          ┃  ┃  ┃    ┃    ┏━━━┻━━━┓     ┃
//     //     31          ┃  ┃  ┃    ┃    ┃       ┃     ┃
//     //   ┏━━┻━━┓       ┃  ┃  ┃    ┃    ┃       ┃     ┃
//     //   ┃     ┃       ┃  ┃ 30    ┃    ┃       ┃     ┃
//     //   ┃     ┃       ┃  ┃ ┏┻━┓  ┃    ┃       ┃     ┃
//     //  29     ┃       ┃  ┃ ┃  ┃  ┃    ┃       ┃     ┃
//     // ┏━┻┓    ┃       ┃  ┃ ┃  ┃  ┃    ┃       ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃   28       ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┏━━┻━┓     ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃   27     ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃  ┏━┻━┓   ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ 26   ┃   ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┏┻━┓ ┃   ┃     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  25     ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┏┻━┓   ┃
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  24
//     // ┃  ┃    ┃       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┏┻━┓
//     // ┃  ┃   23       ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
//     // ┃  ┃ ┏━━┻━┓     ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
//     // ┃  ┃ ┃   22     ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
//     // ┃  ┃ ┃  ┏━┻━━┓  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃  ┃  ┃  ┃  ┃
//     // ┃  ┃ ┃  ┃    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ 21  ┃  ┃  ┃
//     // ┃  ┃ ┃  ┃    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┏┻┓ ┃  ┃  ┃
//     // ┃  ┃ ┃ 20    ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┃ ┃ ┃  ┃  ┃
//     // ┃  ┃ ┃ ┏┻━┓  ┃  ┃  ┃ ┃  ┃  ┃ ┃ ┃  ┃ ┃ ┃ ┃ ┃  ┃  ┃
//     // 0 11 4 7 13 17 12 10 1 18 16 2 6 14 9 3 8 5 15 19

//     std::string const ts_file = "data/allele-frequency-spectrum-simple-example-2.trees";
//     TSKitTreeSequence tree_sequence(ts_file);
//     ForestCompressor  forest_compressor(tree_sequence);
//     CompressedForest  forest = forest_compressor.compress();

//     forest.compute_num_samples_below();
//     // forest.compute_samples_below();

//     std::vector<tsk_id_t> samples     = {0, 11, 4, 7, 13, 17, 12, 10, 1, 18, 16, 2, 6, 14, 9, 3, 8, 5, 15, 19};
//     std::vector<NodeId>   dag_samples = forest_compressor.ts_node2cf_subtree(0, samples);
//     CHECK(dag_samples.size() == samples.size());

//     CHECK(forest.num_nodes() == 39);
//     CHECK(forest.num_trees() == 1);
//     CHECK(forest.num_samples() == samples.size());

//     CHECK(forest.postorder_edges().num_nodes() == forest.num_nodes());
//     CHECK(forest.postorder_edges().num_trees() == forest.num_trees());

//     // The mapping from tree sequence nodes to forest nodes is unique if there's only one tree.
//     std::vector<tsk_id_t> all_nodes(forest.num_nodes());
//     std::iota(all_nodes.begin(), all_nodes.end(), 0);
//     std::vector<NodeId> dag_all_nodes = forest_compressor.ts_node2cf_subtree(0, all_nodes);
//     CHECK(dag_all_nodes.size() == all_nodes.size());

//     auto const num_nodes = dag_all_nodes.size();
//     std::sort(dag_all_nodes.begin(), dag_all_nodes.end());
//     dag_all_nodes.erase(std::unique(dag_all_nodes.begin(), dag_all_nodes.end()), dag_all_nodes.end());
//     CHECK(dag_all_nodes.size() == num_nodes);

//     for (tsk_id_t ts_node = 0; ts_node <= 38; ++ts_node) {
//         auto const dag_node = forest_compressor.ts_node2cf_subtree(0, ts_node);
//         CHECK(dag_node <= 38);
//     }
// }

// TEST_CASE("CompressedForest::ts_nodes2cf_subtree() III", "[CompressedForest]") {
//     std::string const ts_file = "data/allele-frequency-spectrum-simple-example-6.trees";
//     TSKitTreeSequence tree_sequence(ts_file);
//     ForestCompressor  forest_compressor(tree_sequence);
//     CompressedForest  forest = forest_compressor.compress();
//     forest.compute_num_samples_below();
//     // forest.compute_samples_below();
//     std::vector<tsk_id_t> samples     = {0, 8, 11, 6, 5, 19, 12, 13, 15, 1, 2, 18, 4, 7, 14, 16, 17, 9, 10, 3};
//     std::vector<NodeId>   dag_samples = forest_compressor.ts_node2cf_subtree(0, samples);

//     // The samples of both trees should be mapped to the same nodes in the compressed forest DAG.
//     CHECK_THAT(dag_samples, UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(1, samples)));

//     // Each tree has its own root in the compressed forest, even if they would be exactly identical (they aren't in this
//     // example). This simplifies writing code using the compressed forest as the mapping tree <-> root in DAG is
//     // bijective.
//     tsk_id_t const ts_root = 39;
//     CHECK(forest_compressor.ts_node2cf_subtree(0, ts_root) != forest_compressor.ts_node2cf_subtree(1, ts_root));
//     CHECK(dag_samples.size() == samples.size());

//     CHECK(forest.num_trees() == 2);
//     CHECK(forest.num_samples() == samples.size());

//     CHECK(forest.postorder_edges().num_nodes() == forest.num_nodes());
//     CHECK(forest.postorder_edges().num_trees() == forest.num_trees());

//     std::vector<tsk_id_t> all_nodes(tree_sequence.num_nodes());
//     std::iota(all_nodes.begin(), all_nodes.end(), 0);
//     std::vector<NodeId> dag_all_nodes = forest_compressor.ts_node2cf_subtree(0, all_nodes);
//     CHECK(dag_all_nodes.size() == all_nodes.size());

//     auto const num_nodes = dag_all_nodes.size();
//     std::sort(dag_all_nodes.begin(), dag_all_nodes.end());
//     dag_all_nodes.erase(std::unique(dag_all_nodes.begin(), dag_all_nodes.end()), dag_all_nodes.end());
//     CHECK(dag_all_nodes.size() <= num_nodes);

//     for (tsk_id_t ts_node = 0; ts_node <= 39 * 2; ++ts_node) {
//         auto const dag_node = forest_compressor.ts_node2cf_subtree(0, ts_node);
//         CHECK(dag_node < 39 * 2);
//     }
// }

// TODO Think about how we can rewrite these test-cases to not rely on the mapping to be stored
// TEST_CASE("CompressedForest::{num_,}samples_below() IV", "[CompressedForest]") {
//     tsk_treeseq_t tskit_tree_sequence;

//     tsk_treeseq_from_text(
//         &tskit_tree_sequence,
//         10,
//         multi_tree_back_recurrent_nodes,
//         multi_tree_back_recurrent_edges,
//         NULL,
//         multi_tree_back_recurrent_sites,
//         multi_tree_back_recurrent_mutations,
//         multi_tree_back_recurrent_individuals,
//         NULL,
//         0
//     );

//     TSKitTreeSequence tdt_tree_sequence(tskit_tree_sequence); // Takes ownership
//     ForestCompressor  forest_compressor(tdt_tree_sequence);
//     CompressedForest  forest = forest_compressor.compress();
//     forest.compute_num_samples_below();
//     // forest.compute_samples_below();

//     std::vector<tsk_id_t> samples = {0, 1, 2, 3};
//     TreeId const          tree0   = 0;
//     TreeId const          tree1   = 1;
//     TreeId const          tree2   = 2;

//     std::vector<NodeId> dag_samples0 = forest_compressor.ts_node2cf_subtree(tree0, samples);
//     std::vector<NodeId> dag_samples1 = forest_compressor.ts_node2cf_subtree(tree1, samples);
//     std::vector<NodeId> dag_samples2 = forest_compressor.ts_node2cf_subtree(tree2, samples);

//     // The samples of all three trees should be mapped to the same nodes in the compressed forest DAG.
//     CHECK_THAT(dag_samples0, UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree1, samples)));
//     CHECK_THAT(dag_samples1, UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree2, samples)));
//     CHECK_THAT(dag_samples2, UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree0, samples)));

//     // Each tree has its own root in the compressed forest, even if they would be exactly identical (they aren't in this
//     // example). This simplifies writing code using the compressed forest as the mapping tree <-> root in DAG is
//     // bijective.
//     tsk_id_t const ts_root0 = 8;
//     tsk_id_t const ts_root1 = 6;
//     tsk_id_t const ts_root2 = 7;
//     CHECK(forest_compressor.ts_node2cf_subtree(tree0, ts_root0) != forest_compressor.ts_node2cf_subtree(tree1, ts_root1));
//     CHECK(forest_compressor.ts_node2cf_subtree(tree1, ts_root1) != forest_compressor.ts_node2cf_subtree(tree2, ts_root2));
//     CHECK(forest_compressor.ts_node2cf_subtree(tree1, ts_root1) != forest_compressor.ts_node2cf_subtree(tree0, ts_root0));

//     // NodeId const cf_root0 = forest_compressor.ts_node2cf_subtree(tree0, ts_root0);
//     // NodeId const cf_root1 = forest_compressor.ts_node2cf_subtree(tree1, ts_root1);
//     // NodeId const cf_root2 = forest_compressor.ts_node2cf_subtree(tree2, ts_root2);
//     // CHECK_THAT(forest.samples_below(cf_root0), UnorderedRangeEquals(dag_samples0));
//     // CHECK_THAT(forest.samples_below(cf_root1), UnorderedRangeEquals(dag_samples1));
//     // CHECK_THAT(forest.samples_below(cf_root2), UnorderedRangeEquals(dag_samples2));

//     for (NodeId dag_sample: dag_samples0) {
//         CHECK(forest.is_sample(dag_sample));
//         // CHECK_THAT(forest.samples_below(dag_sample), UnorderedRangeEquals(std::vector<NodeId>{dag_sample}));
//     }

//     // CHECK_THAT(
//     //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree0, 6)),
//     //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(0, std::vector<NodeId>{0, 1, 3}))
//     // );
//     // CHECK_THAT(
//     //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree1, 5)),
//     //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree1, std::vector<NodeId>{1, 2, 3}))
//     // );
//     // CHECK_THAT(
//     //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree1, 4)),
//     //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree1, std::vector<NodeId>{2, 3}))
//     // );
//     // CHECK_THAT(
//     //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree2, 5)),
//     //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree2, std::vector<NodeId>{1, 2, 3}))
//     // );
//     // CHECK_THAT(
//     //     forest.samples_below(forest_compressor.ts_node2cf_subtree(tree2, 4)),
//     //     UnorderedRangeEquals(forest_compressor.ts_node2cf_subtree(tree2, std::vector<NodeId>{2, 3}))
//     // );

//     // num_samples_below() for tree 0
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree0, 0)) == 1);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree0, 1)) == 1);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree0, 2)) == 1);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree0, 3)) == 1);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree0, 5)) == 2);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree0, 6)) == 3);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree0, 8)) == 4);

//     // num_samples_below() for tree 1
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree1, 0)) == 1);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree1, 1)) == 1);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree1, 2)) == 1);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree1, 3)) == 1);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree1, 4)) == 2);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree1, 5)) == 3);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree1, 6)) == 4);

//     // num_samples_below() for tree 2
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree2, 0)) == 1);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree2, 1)) == 1);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree2, 2)) == 1);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree2, 3)) == 1);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree2, 4)) == 2);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree2, 5)) == 3);
//     CHECK(forest.num_samples_below(forest_compressor.ts_node2cf_subtree(tree2, 7)) == 4);
// }
