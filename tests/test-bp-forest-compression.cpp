#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_quantifiers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>

#include "mocks/TsToSfMappingExtractor.hpp"
#include "sfkit/SuccinctForest.hpp"
#include "sfkit/assertion_levels.hpp"
#include "sfkit/graph/AdjacencyArrayGraph.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/graph/common.hpp"
#include "sfkit/load/BPForestCompressor.hpp"
#include "sfkit/load/ForestCompressor.hpp"
#include "sfkit/utils/concepts.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace Catch::Matchers;

TEST_CASE("BP Forest Compression Example I", "[BPForestCompression]") {
    /*          6          */
    /*         / \         */
    /*        /   \        */
    /*       /     \       */
    /*      /       5      */
    /*     4       / \     */
    /*    / \     /   \    */
    /*   0   1   2     3   */
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

    TSKitTreeSequence     tree_sequence(tskit_tree_sequence); // Takes ownership
    Ts2SfMappingExtractor ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    auto const            forest = BPForestCompressor(tree_sequence).compress(ts_2_sf_node);

    // CHECK(forest.num_trees() == 1);
    CHECK(forest.num_samples() == 4);
    CHECK(forest.num_leaves() == forest.num_samples());
    // TODO Re-add
    // CHECK(forest.num_unique_subtrees() == 7);
    // CHECK(forest.num_nodes() == forest.num_unique_subtrees());

    CHECK_THAT(forest.is_reference(), NoneTrue());
    CHECK_THAT(forest.balanced_parenthesis(), RangeEquals(std::vector<bool>{1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0,
    0})); CHECK(forest.references().size() == 0);

    // Do not free the tskit tree sequence, as we transferred the ownershop to sfkit_tree_sequence.
    // tsk_treeseq_free(&tskit_tree_sequence);
}

TEST_CASE("BP Forest Compression Example II", "[BPForestCompression]") {
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

    TSKitTreeSequence     tree_sequence(tskit_tree_sequence); // Takes ownership
    Ts2SfMappingExtractor ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    auto const            forest = BPForestCompressor(tree_sequence).compress(ts_2_sf_node);

    // TODO Re-add
    // CHECK(forest.num_trees() == 3);
    CHECK(forest.num_samples() == 4);
    // Even though the last two trees are identical, the root exists twice (see forest compression).
    // TODO Re-add
    // CHECK(forest.num_unique_subtrees() == 11);

    CHECK_THAT(
        forest.is_reference(),
        RangeEquals(std::vector<bool>{
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // tree 0
            0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, // tree 1
            1, 1                                      // tree 2
        })
    );
    CHECK_THAT(
        forest.balanced_parenthesis(),
        RangeEquals(std::vector<bool>{
            1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, // tree 0
            1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, // tree 1
            1, 0                                      // tree 2
        })
    );

    // Each of the four samples is referenced in the second tree. The whole second tree is referenced in the third tree.
    // Each reference consists of two elements: <start, length>.
    CHECK(forest.references().size() == 10);
}

// test-zazu.trees
// 76.16┊         ┊         ┊         ┊         ┊   9     ┊ 
//      ┊         ┊         ┊         ┊         ┊ ┏━┻━┓   ┊ 
// 50.34┊   8     ┊         ┊   8     ┊    8    ┊ ┃   8   ┊ 
//      ┊ ┏━┻━┓   ┊         ┊ ┏━┻━┓   ┊  ┏━┻━┓  ┊ ┃  ┏┻━┓ ┊ 
// 27.71┊ ┃   ┃   ┊         ┊ ┃   ┃   ┊  7   ┃  ┊ ┃  ┃  ┃ ┊ 
//      ┊ ┃   ┃   ┊         ┊ ┃   ┃   ┊ ┏┻┓  ┃  ┊ ┃  ┃  ┃ ┊ 
// 21.39┊ ┃   ┃   ┊   6     ┊ ┃   ┃   ┊ ┃ ┃  ┃  ┊ ┃  ┃  ┃ ┊ 
//      ┊ ┃   ┃   ┊ ┏━┻━┓   ┊ ┃   ┃   ┊ ┃ ┃  ┃  ┊ ┃  ┃  ┃ ┊ 
// 11.54┊ ┃   5   ┊ ┃   5   ┊ ┃   5   ┊ ┃ ┃  5  ┊ ┃  5  ┃ ┊ 
//      ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┏┻┓ ┃ ┊ 
// 5.30 ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ 
//      ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ 
// 0.00 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊ 0 3 1 2 ┊ 0 1 2 3 ┊ 
//      0         2        11        13        19        20 

// test-pumbaa.trees
// 107.57┊    15   ┊   15    ┊  15     ┊  15     ┊  15     ┊  15     ┊  15     ┊         ┊         ┊         ┊         ┊ 
//       ┊   ┏━┻━┓ ┊  ┏━┻━┓  ┊ ┏━┻━┓   ┊ ┏━┻━┓   ┊ ┏━┻━┓   ┊ ┏━┻━┓   ┊ ┏━┻━┓   ┊         ┊         ┊         ┊         ┊ 
// 54.45 ┊  14   ┃ ┊  ┃   ┃  ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊         ┊         ┊         ┊         ┊ 
//       ┊  ┏┻━┓ ┃ ┊  ┃   ┃  ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊         ┊         ┊         ┊         ┊ 
// 32.13 ┊  ┃  ┃ ┃ ┊  ┃   ┃  ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊         ┊         ┊    13   ┊   13    ┊ 
//       ┊  ┃  ┃ ┃ ┊  ┃   ┃  ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊         ┊         ┊   ┏━┻━┓ ┊  ┏━┻━┓  ┊ 
// 29.26 ┊  ┃  ┃ ┃ ┊  ┃  12  ┊ ┃  12   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊         ┊         ┊   ┃   ┃ ┊  ┃   ┃  ┊ 
//       ┊  ┃  ┃ ┃ ┊  ┃  ┏┻┓ ┊ ┃  ┏┻━┓ ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊         ┊         ┊   ┃   ┃ ┊  ┃   ┃  ┊ 
// 29.02 ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊ ┃  ┃  ┃ ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃  11   ┊    11   ┊    11   ┊   ┃   ┃ ┊  ┃   ┃  ┊ 
//       ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊ ┃  ┃  ┃ ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃   ┃   ┊ ┃ ┏━┻┓  ┊   ┏━┻━┓ ┊   ┏━┻━┓ ┊   ┃   ┃ ┊  ┃   ┃  ┊ 
// 29.01 ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊ ┃ 10  ┃ ┊ ┃  10   ┊ ┃  10   ┊ ┃   ┃   ┊ ┃ ┃  ┃  ┊   ┃   ┃ ┊   ┃   ┃ ┊   ┃   ┃ ┊  ┃   ┃  ┊ 
//       ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊ ┃ ┏┻┓ ┃ ┊ ┃  ┏┻━┓ ┊ ┃  ┏┻━┓ ┊ ┃   ┃   ┊ ┃ ┃  ┃  ┊   ┃   ┃ ┊   ┃   ┃ ┊   ┃   ┃ ┊  ┃   ┃  ┊ 
// 23.45 ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃  9  ┃ ┊ ┃  ┃  ┃ ┊ ┃   ┃   ┊ ┃ ┃  ┃  ┊   9   ┃ ┊   9   ┃ ┊   9   ┃ ┊  ┃   9  ┊ 
//       ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┏┻┓ ┃ ┊ ┃  ┃  ┃ ┊ ┃   ┃   ┊ ┃ ┃  ┃  ┊ ┏━┻┓  ┃ ┊  ┏┻━┓ ┃ ┊  ┏┻━┓ ┃ ┊  ┃  ┏┻┓ ┊ 
// 22.08 ┊  8  ┃ ┃ ┊  8  ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃  ┃  ┃ ┊ ┃   ┃   ┊ ┃ ┃  ┃  ┊ ┃  ┃  ┃ ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊ 
//       ┊ ┏┻┓ ┃ ┃ ┊ ┏┻┓ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃  ┃  ┃ ┊ ┃   ┃   ┊ ┃ ┃  ┃  ┊ ┃  ┃  ┃ ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊ 
// 14.42 ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃  7  ┃ ┊ ┃   7   ┊ ┃ ┃  ┃  ┊ ┃  ┃  ┃ ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊ 
//       ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┏┻┓ ┃ ┊ ┃ ┏━┻┓  ┊ ┃ ┃  ┃  ┊ ┃  ┃  ┃ ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊ 
// 12.42 ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃  ┃  ┊ ┃ ┃  ┃  ┊ ┃  ┃  ┃ ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊  6  ┃ ┃ ┊ 
//       ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃  ┃  ┊ ┃ ┃  ┃  ┊ ┃  ┃  ┃ ┊  ┃  ┃ ┃ ┊  ┃  ┃ ┃ ┊ ┏┻┓ ┃ ┃ ┊ 
// 9.06  ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃  ┃  ┊ ┃ ┃  ┃  ┊ ┃  ┃  ┃ ┊  5  ┃ ┃ ┊  5  ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ 
//       ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃  ┃  ┊ ┃ ┃  ┃  ┊ ┃  ┃  ┃ ┊ ┏┻┓ ┃ ┃ ┊ ┏┻┓ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ 
// 5.32  ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊ ┃  4  ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ 
//       ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┏┻┓ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊ 
// 0.00  ┊ 0 2 1 3 ┊ 0 2 1 3 ┊ 0 1 2 3 ┊ 0 1 3 2 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊ 0 2 3 1 ┊ 0 3 2 1 ┊ 0 3 2 1 ┊ 0 1 2 3 ┊ 
//       0        10        11        17        18        23        32        38        49        68        93        100
