#include "tskit/trees.h"
#include <debug/vector>
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <stddef.h>

#include "catch2/catch_assertion_info.hpp"
#include "sfkit/graph/primitives.hpp"
#include "sfkit/tskit/tskit.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch;
using namespace ::Catch::Matchers;

using namespace sfkit::graph;
using namespace sfkit::tskit;

TEST_CASE("TSKitTree::eulertour() Example I", "[TSKitTree]") {
    //     8   ┊         ┊         ┊
    //   ┏━┻━┓ ┊         ┊         ┊
    //   ┃   ┃ ┊         ┊   7     ┊
    //   ┃   ┃ ┊         ┊ ┏━┻━┓   ┊
    //   6   ┃ ┊   6     ┊ ┃   ┃   ┊
    // ┏━┻┓  ┃ ┊ ┏━┻━┓   ┊ ┃   ┃   ┊
    // ┃  5  ┃ ┊ ┃   5   ┊ ┃   5   ┊
    // ┃ ┏┻┓ ┃ ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊
    // ┃ ┃ ┃ ┃ ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊
    // ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊
    // 0 1 3 2 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊

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

    TSKitTree tree{sfkit_tree_sequence};

    tree.first();
    CHECK_THAT(tree.eulertour(), RangeEquals(std::vector<NodeId>{8, 2, 6, 0, 5, 1, 3, 5, 6, 8}));

    tree.next();
    CHECK_THAT(tree.eulertour(), RangeEquals(std::vector<NodeId>{6, 0, 5, 1, 4, 2, 3, 4, 5, 6}));

    tree.next();
    CHECK_THAT(tree.eulertour(), RangeEquals(std::vector<NodeId>{7, 0, 5, 1, 4, 2, 3, 4, 5, 7}));
}

TEST_CASE("TSKitTree::eulertour() Example II", "[TSKitTree]") {
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

    TSKitTreeSequence sfkit_tree_sequence(std::move(tskit_tree_sequence));
    REQUIRE(sfkit_tree_sequence.is_owning());

    TSKitTree tree{sfkit_tree_sequence};

    tree.first();
    CHECK_THAT(tree.eulertour(), RangeEquals(std::vector<NodeId>{6, 4, 0, 1, 4, 5, 2, 3, 5, 6}));
}

TEST_CASE("TSKitTree::eulertour() Timon", "[TSKitTree]") {
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

    std::string const ts_file = "data/test-timon.trees";
    TSKitTreeSequence tree_sequence(ts_file);

    TSKitTree tree{tree_sequence};

    tree.first();
    CHECK_THAT(tree.eulertour(), RangeEquals(std::vector<NodeId>{38, 34, 17, 30, 15, 16, 30, 34, 37, 35, 3,  5,
                                                                 35, 36, 22, 18, 20, 1,  2,  20, 22, 33, 27, 19,
                                                                 26, 4,  21, 7,  14, 21, 26, 27, 32, 29, 8,  28,
                                                                 13, 25, 12, 23, 9,  11, 23, 25, 28, 29, 31, 10,
                                                                 24, 0,  6,  24, 31, 32, 33, 36, 37, 38}));
}

TEST_CASE("TSKitTree::eulertour() Scar", "[TSKitTree]") {
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

    std::string const ts_file = "data/test-scar.trees";
    TSKitTreeSequence tree_sequence(ts_file);

    TSKitTree tree{tree_sequence};

    tree.first();
    CHECK_THAT(tree.eulertour(), RangeEquals(std::vector<NodeId>{38, 33, 24, 15, 19, 24, 32, 25, 5,  21, 3,  8,
                                                                 21, 25, 28, 2,  27, 9,  26, 6,  14, 26, 27, 28,
                                                                 32, 33, 37, 35, 16, 30, 1,  18, 30, 35, 36, 10,
                                                                 34, 12, 31, 23, 4,  22, 17, 20, 7,  13, 20, 22,
                                                                 23, 29, 0,  11, 29, 31, 34, 36, 37, 38}));
}

TEST_CASE("TSKitTree::eulertour() Shenzi", "[TSKitTree]") {
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

    std::string const ts_file = "data/test-shenzi.trees";
    TSKitTreeSequence tree_sequence(ts_file);

    TSKitTree tree{tree_sequence};

    tree.first();
    CHECK_THAT(tree.eulertour(), RangeEquals(std::vector<NodeId>{38, 33, 24, 15, 19, 24, 32, 25, 5,  21, 3,  8,
                                                                 21, 25, 28, 2,  27, 9,  26, 6,  14, 26, 27, 28,
                                                                 32, 33, 37, 35, 16, 30, 1,  18, 30, 35, 36, 10,
                                                                 34, 12, 31, 23, 4,  22, 17, 20, 7,  13, 20, 22,
                                                                 23, 29, 0,  11, 29, 31, 34, 36, 37, 38}));
    // Do not free the tskit tree sequence, as we transferred the ownershop to sfkit_tree_sequence.
    // tsk_treeseq_free(&tskit_tree_sequence);
}

TEST_CASE("TSKitTree::invalidated_nodes() Example I", "[TSKitTree]") {
    //     8   ┊         ┊         ┊
    //   ┏━┻━┓ ┊         ┊         ┊
    //   ┃   ┃ ┊         ┊   7     ┊
    //   ┃   ┃ ┊         ┊ ┏━┻━┓   ┊
    //   6   ┃ ┊   6     ┊ ┃   ┃   ┊
    // ┏━┻┓  ┃ ┊ ┏━┻━┓   ┊ ┃   ┃   ┊
    // ┃  5  ┃ ┊ ┃   5   ┊ ┃   5   ┊
    // ┃ ┏┻┓ ┃ ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊
    // ┃ ┃ ┃ ┃ ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊
    // ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊
    // 0 1 3 2 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊

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

    TSKitTree tree{sfkit_tree_sequence};

    tree.first();

    tree.next();
    CHECK_THAT(tree.invalidated_nodes(), UnorderedRangeEquals(std::vector<NodeId>{4, 5, 6, 8}));

    tree.next();
    CHECK_THAT(tree.invalidated_nodes(), UnorderedRangeEquals(std::vector<NodeId>{6, 7}));
}

TEST_CASE("TSKitTree::invalidated_nodes() Zazu", "[BPForestCompression]") {
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

    TSKitTreeSequence tree_sequence("data/test-zazu.trees");
    REQUIRE(tree_sequence.is_owning());

    TSKitTree tree{tree_sequence};

    tree.first();

    // Tree 0 -> 1
    tree.next();
    CHECK_THAT(tree.invalidated_nodes(), UnorderedRangeEquals(std::vector<NodeId>{6, 8}));

    // Tree 1 -> 2
    tree.next();
    CHECK_THAT(tree.invalidated_nodes(), UnorderedRangeEquals(std::vector<NodeId>{6, 8}));

    // Tree 2 -> 3
    tree.next();
    CHECK_THAT(tree.invalidated_nodes(), UnorderedRangeEquals(std::vector<NodeId>{4, 5, 7, 8}));

    // Tree 3 -> 4
    tree.next();
    CHECK_THAT(tree.invalidated_nodes(), UnorderedRangeEquals(std::vector<NodeId>{7, 8, 9}));
}
