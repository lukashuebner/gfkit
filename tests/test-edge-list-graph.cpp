#include <algorithm>
#include <debug/vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <debug/safe_iterator.h>
#include <stddef.h>

#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/graph/primitives.hpp"

using namespace Catch::Matchers;

using namespace sfkit::graph;

TEST_CASE("EdgeListGraph Basics", "[EdgeListGraph]") {
    EdgeListGraph graph;

    SECTION("Empty graph") {
        CHECK(graph.num_edges() == 0);
        CHECK(graph.begin() == graph.end());
    }

    SECTION("Graph with a single edge") {
        graph.insert_edge(0, 1);
        graph.compute_num_nodes();
        CHECK(graph.num_nodes() == 2);
        CHECK(graph.num_edges() == 1);
        CHECK(graph.begin() != graph.end());
        CHECK(graph.begin()->from() == 0);
        CHECK(graph.begin()->from() == 0); // Calling this twice should not change the result
        CHECK(graph.begin()->to() == 1);
        CHECK(graph.begin()->to() == 1); // Calling this twice should not change the result
    }

    SECTION("Set the number of nodes to a wrong value") {
        // If we add a single edge, the number of nodes should be autocomputed in O(|edges|) time.
        graph.insert_edge(0, 1);
        graph.compute_num_nodes();
        CHECK(graph.num_nodes() == 2);
        CHECK(graph.num_edges() == 1);

        // Invalidate the num_nodes cache and check that the value is recomputed.
        graph.insert_edge(0, 2);
        graph.compute_num_nodes();
        CHECK(graph.num_nodes() == 3);
    }

    SECTION("Graph with three edges") {
        graph.insert_edge(0, 2);
        graph.insert_edge(1, 2);
        graph.insert_edge(2, 2);
        graph.compute_num_nodes();
        CHECK(graph.num_edges() == 3);
    }

    SECTION("Graph with multiedges") {
        // We do not check that an edge is unique (too slow), so adding the same edge twice (or three times) should
        // work.
        graph.insert_edge(0, 1);
        graph.insert_edge(0, 1);
        graph.compute_num_nodes();
        CHECK(graph.num_edges() == 2);
    }
}

TEST_CASE("EdgeListGraph::add_leaves()", "[EdgeListGraph]") {
    EdgeListGraph       graph;
    std::vector<NodeId> leaves = {};

    CHECK(graph.leaves().size() == 0);
    CHECK(graph.leaves() == leaves);

    graph.insert_leaf(0);
    leaves.push_back(0);
    CHECK(graph.leaves().size() == 1);
    CHECK(graph.leaves() == leaves);

    // No interference between leaves and roots.
    graph.insert_root(12);
    CHECK(graph.leaves().size() == 1);
    CHECK(graph.leaves() == leaves);

    for (NodeId i: std::vector<NodeId>{1, 2, 5}) {
        graph.insert_leaf(i);
        leaves.push_back(i);
    }
    CHECK(graph.leaves().size() == 4);
    CHECK_THAT(graph.leaves(), RangeEquals(leaves));
}

TEST_CASE("EdgeListGraph::add_root()", "[EdgeListGraph]") {
    EdgeListGraph       graph;
    std::vector<NodeId> roots = {};

    CHECK(graph.roots().size() == 0);
    CHECK(graph.roots() == roots);

    graph.insert_root(0);
    roots.push_back(0);
    CHECK(graph.roots().size() == 1);
    CHECK(graph.roots() == roots);

    // No interference between leaves and roots.
    graph.insert_leaf(12);
    CHECK(graph.roots().size() == 1);
    CHECK(graph.roots() == roots);

    for (NodeId i: std::vector<NodeId>{1, 2, 5}) {
        graph.insert_root(i);
        roots.push_back(i);
    }
    CHECK(graph.roots().size() == 4);
    CHECK_THAT(graph.roots(), RangeEquals(roots));
}

TEST_CASE("EdgeListGraph::directed()", "[EdgeListGraph]") {
    EdgeListGraph graph;

    CHECK(graph.directed() == true);

    graph.insert_edge(0, 1);
    CHECK(graph.directed() == true);

    graph.insert_edge(1, 0);
    CHECK(graph.directed() == true);
}

TEST_CASE("EdgeListGraph::traversal_order()", "[EdgeListGraph]") {
    EdgeListGraph graph;

    CHECK(graph.traversal_order() == TraversalOrder::Unordered);

    auto traversal_orders = std::vector<TraversalOrder>{
        TraversalOrder::Unordered,
        TraversalOrder::Inorder,
        TraversalOrder::Preorder,
        TraversalOrder::Postorder,
        TraversalOrder::Levelorder};

    for (TraversalOrder order: traversal_orders) {
        graph.traversal_order(order);
        CHECK(graph.traversal_order() == order);
        switch (order) {
            case TraversalOrder::Unordered:
                CHECK(graph.is_unordered());
                CHECK_FALSE(graph.is_preorder());
                CHECK_FALSE(graph.is_inorder());
                CHECK_FALSE(graph.is_postorder());
                CHECK_FALSE(graph.is_levelorder());
                break;
            case TraversalOrder::Preorder:
                CHECK_FALSE(graph.is_unordered());
                CHECK(graph.is_preorder());
                CHECK_FALSE(graph.is_inorder());
                CHECK_FALSE(graph.is_postorder());
                CHECK_FALSE(graph.is_levelorder());
                break;
            case TraversalOrder::Inorder:
                CHECK_FALSE(graph.is_unordered());
                CHECK_FALSE(graph.is_preorder());
                CHECK(graph.is_inorder());
                CHECK_FALSE(graph.is_postorder());
                CHECK_FALSE(graph.is_levelorder());
                break;
            case TraversalOrder::Postorder:
                CHECK_FALSE(graph.is_unordered());
                CHECK_FALSE(graph.is_preorder());
                CHECK_FALSE(graph.is_inorder());
                CHECK(graph.is_postorder());
                CHECK_FALSE(graph.is_levelorder());
                break;
            case TraversalOrder::Levelorder:
                CHECK_FALSE(graph.is_unordered());
                CHECK_FALSE(graph.is_preorder());
                CHECK_FALSE(graph.is_inorder());
                CHECK_FALSE(graph.is_postorder());
                CHECK(graph.is_levelorder());
                break;
        }
    }
}

TEST_CASE("EdgeListGraph::edges_are_sorted()", "[EdgeListGraph]") {
    using SortBy = EdgeListGraph::SortBy;

    EdgeListGraph graph;

    SECTION("Empty graph") {
        CHECK(graph.edges_are_sorted(SortBy::FromVertex));
        CHECK(graph.edges_are_sorted(SortBy::ToVertex));
    }

    SECTION("Consecutive, both sorted") {
        for (NodeId i: std::vector<NodeId>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}) {
            graph.insert_edge(i, i + 1);
        }
        CHECK(graph.edges_are_sorted(SortBy::FromVertex));
        CHECK(graph.edges_are_sorted(SortBy::ToVertex));
    }

    SECTION("Non-consecutive, both sorted") {
        for (NodeId i: std::vector<NodeId>{0, 10, 21, 34, 40, 544}) {
            graph.insert_edge(i, i + 1);
        }
        CHECK(graph.edges_are_sorted(SortBy::FromVertex));
        CHECK(graph.edges_are_sorted(SortBy::ToVertex));
    }

    SECTION("Both unsorted") {
        for (NodeId i: std::vector<NodeId>{0, 3, 2, 30, 4, 5, 6, 7, 8, 9}) {
            graph.insert_edge(i, i + 1);
        }
        CHECK_FALSE(graph.edges_are_sorted(SortBy::FromVertex));
        CHECK_FALSE(graph.edges_are_sorted(SortBy::ToVertex));
    }

    SECTION("Both sorted with double edges") {
        for (NodeId i: std::vector<NodeId>{0, 10, 21, 34, 40, 544}) {
            graph.insert_edge(i, i + 1);
            graph.insert_edge(i, i + 1);
        }
        CHECK(graph.edges_are_sorted(SortBy::FromVertex));
        CHECK(graph.edges_are_sorted(SortBy::ToVertex));
    }

    SECTION("sorted by from, not by to vertex") {
        graph.insert_edge(0, 1);
        graph.insert_edge(1, 4);
        graph.insert_edge(2, 1);
        graph.insert_edge(5, 0);
        graph.insert_edge(8, 1);
    }

    SECTION("sorted by to, not by from vertex") {
        graph.insert_edge(10, 1);
        graph.insert_edge(1, 4);
        graph.insert_edge(2, 6);
        graph.insert_edge(34, 6);
        graph.insert_edge(8, 8);
    }
}

TEST_CASE("EdgeListGraph::sort() simple graph", "[EdgeListGraph]") {
    using SortBy = EdgeListGraph::SortBy;

    EdgeListGraph graph;
    size_t        num_edges = 0;
    for (NodeId i: std::vector<NodeId>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}) {
        graph.insert_edge(i, i + 1);
        num_edges++;
    }
    for (NodeId i: std::vector<NodeId>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}) {
        graph.insert_edge(i + 2, i);
        num_edges++;
    }
    CHECK(graph.num_edges() == num_edges);

    SECTION("Sort edges by from vertex") {
        CHECK_FALSE(graph.edges_are_sorted(SortBy::FromVertex));
        CHECK_FALSE(graph.edges_are_sorted(SortBy::ToVertex));

        graph.sort_edges(SortBy::FromVertex);

        CHECK(graph.edges_are_sorted(SortBy::FromVertex));
        CHECK_FALSE(graph.edges_are_sorted(SortBy::ToVertex));
    }

    SECTION("Sort edges by to vertex") {
        CHECK_FALSE(graph.edges_are_sorted(SortBy::FromVertex));
        CHECK_FALSE(graph.edges_are_sorted(SortBy::ToVertex));

        graph.sort_edges(SortBy::ToVertex);

        CHECK(graph.edges_are_sorted(SortBy::ToVertex));
        CHECK_FALSE(graph.edges_are_sorted(SortBy::FromVertex));
    }
}
