
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <numeric>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_quantifiers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <tdt/graph/adjacency-array-graph.hpp>
#include <tdt/graph/common.hpp>
#include <tdt/graph/edge-list-graph.hpp>
#include <tskit.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/load/forest-compressor.hpp"
#include "tdt/tskit.hpp"

using namespace ::Catch::Matchers;

std::string const ts_file1          = "data/generated-pop100-seq100000000-ind10-seed0001.trees";
std::string const ts_file3          = "data/generated-pop10000-seq10000-ind20-seed0001.trees";
std::string const ts_file4          = "data/generated-pop10000-seq10000000-ind20-seed0001.trees";
std::string const ts_file5          = "data/generated-pop100-seq100000-ind100-seed0001.trees";
std::string const ts_file6          = "data/generated-pop10000-seq100000-ind1000-seed0001.trees";
std::string const tgp_chr8_ts_file  = "data/tgp_chr8.trees";
std::string const tgp_chr14_ts_file = "data/tgp_chr14.trees";
std::string const tgp_chr19_ts_file = "data/tgp_chr19.trees";
std::string const tgp_chr20_ts_file = "data/tgp_chr20.trees";
std::string const tgp_chr21_ts_file = "data/tgp_chr21.trees";
std::string const tgp_chr22_ts_file = "data/tgp_chr22.trees";

// TODO Use custom main to parse the command line arguments and set the expensive tests flag at runtime
#ifdef TDT_EXPENSIVE_TESTS
std::vector<std::string> ts_files = load_file_names_from_file("data/expensive_test_files.txt")
#else
std::vector<std::string> ts_files{ts_file1, ts_file3, ts_file5, ts_file6};
#endif

// TEST_CASE("Build tree sequence and tree from tgp_chr8.trees", "[LoadFromTreeSequence]") {
//     TSKitTreeSequence tree_sequence(tgp_chr8_ts_file);
//     TSKitTree         current_tree(tree_sequence);

//     CHECK(tree_sequence.num_nodes() == 44910);
//     CHECK(tree_sequence.num_edges() == 382351);
//     CHECK(tree_sequence.num_sites() == 106678);
//     CHECK(tree_sequence.num_mutations() == 391346);
//     CHECK(tree_sequence.num_trees() == 9973);
//     CHECK(tree_sequence.num_individuals() == 2504);
//     CHECK(tree_sequence.num_populations() == 26);
//     CHECK(tree_sequence.num_samples() == 5008);
//     CHECK(tree_sequence.num_samples() == tree_sequence.num_individuals() * 2); // humans are diploid
//     CHECK(tree_sequence.sequence_length() == 146364024.0);

//     for (current_tree.first(); current_tree.is_valid(); current_tree.next()) {
//         CHECK(current_tree.num_roots() == 1);
//     }
// }

TEST_CASE("TSKitTree basics", "[LoadFromTreeSequence]") {
    auto const&       ts_file = GENERATE_REF(from_range(ts_files));
    TSKitTreeSequence tree_sequence(ts_file);
    TSKitTree         tree(tree_sequence);

    CHECK(tree.first());
    CHECK(tree.is_valid());
    CHECK(tree.next());
    CHECK(tree.is_valid());

    size_t num_trees = 0;
    for (tree.first(); tree.is_valid(); tree.next()) {
        num_trees++;
        CHECK(tree.is_valid());
        CHECK(tree.is_root(tree.root()));
        CHECK(tree.num_samples() == tree_sequence.num_samples());

        size_t num_roots   = 0;
        size_t num_samples = 0;
        for (auto const node: tree.postorder()) {
            CHECK_FALSE((tree.is_root(node) && tree.is_sample(node)));
            if (tree.is_root(node)) {
                num_roots++;
            } else if (tree.is_sample(node)) {
                num_samples++;
            }

            if (tree.is_sample(node)) {
                CHECK(tree.num_children(node) == 0);
            } else {
                CHECK(tree.num_children(node) > 0);
            }
        }
        CHECK(num_roots == tree.num_roots());
        CHECK(tree.num_roots() == 1);
        CHECK(num_samples == tree_sequence.num_samples());
    }
    CHECK(num_trees == tree_sequence.num_trees());
}

TEST_CASE("TSKitTree.postorder()", "[LoadFromTreeSequence]") {
    auto const&       ts_file = GENERATE_REF(from_range(ts_files));
    TSKitTreeSequence tree_sequence(ts_file);
    TSKitTree         tree(tree_sequence);

    for (tree.first(); tree.is_valid(); tree.next()) {
        auto              postorder = tree.postorder();
        auto const        num_nodes = postorder.size();
        std::vector<bool> visited(tree.max_node_id() + 1, false);
        // visited[asserting_cast<size_t>(tree.root())] = true;
        for (auto const node: postorder) {
            size_t node_idx = asserting_cast<size_t>(node);
            CHECK_FALSE(visited[node_idx]);
            for (auto child: tree.children(node)) {
                size_t child_idx = asserting_cast<size_t>(child);
                CHECK(visited[child_idx]);
            }
            visited[node_idx] = true;
        }
        auto const num_visited = std::accumulate(visited.begin(), visited.end(), 0ul, std::plus<>());
        CHECK(num_visited == num_nodes);
    }
}

TEST_CASE("TSKitTree.children()", "[LoadFromTreeSequence]") {
    auto const&       ts_file = GENERATE_REF(from_range(ts_files));
    TSKitTreeSequence tree_sequence(ts_file);
    TSKitTree         tree(tree_sequence);

    for (tree.first(); tree.is_valid(); tree.next()) {
        std::vector<bool> appeared_as_child(tree.max_node_id() + 1, false);
        appeared_as_child[asserting_cast<size_t>(tree.root())] = true;
        std::vector<bool> postorder_visited(tree.max_node_id() + 1, false);
        auto const        postorder = tree.postorder();
        auto const        num_nodes = postorder.size();
        for (auto const node: postorder) {
            CHECK_FALSE(postorder_visited[asserting_cast<size_t>(node)]);
            postorder_visited[asserting_cast<size_t>(node)] = true;

            auto children = tree.children(node);
            for (auto child: children) {
                CHECK(postorder_visited[asserting_cast<size_t>(child)]);
                CHECK_FALSE(appeared_as_child[asserting_cast<size_t>(child)]);
                appeared_as_child[asserting_cast<size_t>(child)] = true;
                CHECK(tree.parent(child) == node);
            }
        }
        auto const num_appeared_as_child =
            std::accumulate(appeared_as_child.begin(), appeared_as_child.end(), 0ul, std::plus<>());
        auto const num_postorder_visited =
            std::accumulate(postorder_visited.begin(), postorder_visited.end(), 0ul, std::plus<>());
        CHECK(num_appeared_as_child == num_nodes);
        CHECK(num_postorder_visited == num_nodes);
    }
}

// TODO Add sanitizers to test binaries
TEST_CASE("ForestCompressor", "[LoadFromTreeSequence]") {
    auto const& ts_file = GENERATE_REF(from_range(ts_files));

    TSKitTreeSequence tree_sequence(ts_file);
    ForestCompressor  forest_compressor(tree_sequence);
    CompressedForest  compressed_forest = forest_compressor.compress();

    EdgeListGraph const& dag = compressed_forest.postorder_edges();

    REQUIRE(dag.is_postorder());
    CHECK(dag.num_trees() == tree_sequence.num_trees());
    CHECK(dag.num_roots() == tree_sequence.num_trees());
    CHECK(dag.roots().size() == tree_sequence.num_trees());
    CHECK(dag.num_leaves() == tree_sequence.num_samples());

    // Compute the number of samples below each root.
    std::vector<size_t> subtree_sizes(dag.num_nodes(), 0);
    for (NodeId leaf: dag.leaves()) {
        subtree_sizes[leaf] = 1;
    }
    for (auto& edge: dag) {
        subtree_sizes[edge.from()] += subtree_sizes[edge.to()];
    }

    // Check the computed subtee sizes.
    std::for_each(subtree_sizes.begin(), subtree_sizes.end(), [](auto& x) { CHECK(x >= 1); });
    for (NodeId leaf: dag.leaves()) {
        CHECK(subtree_sizes[leaf] == 1);
    }
    for (NodeId root: dag.roots()) {
        CHECK(subtree_sizes[root] == dag.num_leaves());
    }

    // If the edge list is in postorder, then all in edges of a node will appear after all out edges of a node.
    std::vector<bool> pointed_to(dag.num_nodes(), false);
    for (auto& edge: dag) {
        CHECK_FALSE(pointed_to[edge.from()]);
        pointed_to[edge.to()] = true;
    }

    // Do all nodes have in edges or are root nodes?
    for (auto root: dag.roots()) {
        // Roots should not have any in edges.
        CHECK_FALSE(pointed_to[root]);
        pointed_to[root] = true; // This allows us to check if every node has been pointed to more easily.
    }
    auto const num_pointed_to = std::count(pointed_to.begin(), pointed_to.end(), true);
    CHECK(asserting_cast<size_t>(num_pointed_to) == dag.num_nodes());
    CHECK_THAT(pointed_to, AllTrue());
}
