#include <debug/vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include "sfkit/graph/AdjacencyArrayGraph.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/graph/common.hpp"

using namespace Catch::Matchers;

TEST_CASE("AdjacencyArrayGraph Basics", "[AdjacencyArrayGraph]") {
    EdgeListGraph edge_list;
    edge_list.insert_edge(0, 1);
    edge_list.insert_edge(1, 2);
    edge_list.insert_edge(2, 1);
    edge_list.insert_edge(2, 0);
    edge_list.insert_edge(3, 0);
    edge_list.insert_root(3);
    edge_list.compute_num_nodes();

    AdjacencyArrayGraph graph(edge_list);

    CHECK(graph.num_nodes() == 4);
    CHECK(graph.num_edges() == 5);
    CHECK(graph.has_edge(0, 1));
    CHECK(graph.has_edge(1, 2));
    CHECK(graph.has_edge(2, 1));
    CHECK(graph.has_edge(2, 0));
    CHECK(graph.has_edge(3, 0));
    CHECK_FALSE(graph.has_edge(1, 0));
    CHECK_FALSE(graph.has_edge(2, 2));
    CHECK_FALSE(graph.has_edge(0, 0));
    CHECK_FALSE(graph.has_edge(0, 2));
    CHECK_FALSE(graph.has_edge(0, 3));
    CHECK_FALSE(graph.has_edge(1, 3));
    CHECK_FALSE(graph.has_edge(2, 3));
    CHECK_FALSE(graph.has_edge(3, 3));

    CHECK_THAT(graph.adjacent_vertices(0), RangeEquals(std::vector<NodeId>{1}));
    CHECK_THAT(graph.adjacent_vertices(1), RangeEquals(std::vector<NodeId>{2}));
    CHECK_THAT(graph.adjacent_vertices(2), RangeEquals(std::vector<NodeId>{0, 1}));

    CHECK_THAT(graph.roots(), RangeEquals(std::vector<NodeId>{3}));
    CHECK_THAT(graph.leaves(), RangeEquals(std::vector<NodeId>{}));

    // TODO Check that has_edge() and adjacent_vertices() throw assertions for invalid vertex ids.
}
