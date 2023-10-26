#pragma once

#include <vector>

#include <kassert/kassert.hpp>

#include "EdgeListGraph.hpp"
#include "common.hpp"
#include "sfkit/assertion_levels.hpp"
#include "sfkit/checking_casts.hpp"

class AdjacencyArrayGraph {
public:
    // TODO think about moving the roots and leaves
    AdjacencyArrayGraph() = default;
    AdjacencyArrayGraph(EdgeListGraph& edges, bool sorted = false) {
        initialize(edges, sorted);
    }

    // Not using a constructor enables easier serialization
    void initialize(EdgeListGraph& edges, bool sorted = false) {
        _num_edges = 0;
        _roots     = edges.roots();
        _leaves    = edges.leaves();

        if (!sorted) {
            edges.sort_edges(EdgeListGraph::SortBy::FromVertex);
        }

        auto num_vertices = edges.num_nodes();
        _adjacency_array.resize(asserting_cast<uint32_t>(num_vertices));

        for (auto const& edge: edges) {
            _adjacency_array[edge.from()].push_back(edge.to());
            _num_edges++;
        }

        // Sort the adjacency lists such that the to vertices are in ascending order.
        // This enables us to use binary search to check for the existence of an edge.
        // TODO For sparse graphs, maybe a linear search is faster
        std::for_each(_adjacency_array.begin(), _adjacency_array.end(), [](auto& vertex) {
            std::sort(vertex.begin(), vertex.end());
        });

// TODO Check that there are no duplicate edges
#if KASSERT_ENABLED(SFKIT_ASSERTION_LEVEL_NORMAL)
        std::for_each(_adjacency_array.begin(), _adjacency_array.end(), [](auto& vertex) {
            auto last = std::unique(vertex.begin(), vertex.end());
            vertex.erase(last, vertex.end());
        });
#endif

        // TODO Precompute the traversal order for cache efficient traversals
    }

    NodeId num_nodes() const {
        return asserting_cast<NodeId>(_adjacency_array.size());
    }

    NodeId num_edges() const {
        return _num_edges;
    }

    bool has_edge(NodeId from, NodeId to) const {
        KASSERT(from < num_nodes(), "from vertex id greater than the maximum vertex id.", sfkit::assert::light);
        KASSERT(to < num_nodes(), "to vertex id greater than the maximum vertex id.", sfkit::assert::light);

        // TODO For sparse graphs, maybe a linear search is faster
        auto const& vertex = _adjacency_array[from];
        return std::binary_search(vertex.begin(), vertex.end(), to);
    }

    std::vector<NodeId> const& adjacent_vertices(NodeId vertex) const {
        KASSERT(vertex < num_nodes(), "vertex id greater than the maximum vertex id.", sfkit::assert::light);
        return _adjacency_array[vertex];
    }

    std::vector<NodeId> const& roots() const {
        return _roots;
    }

    std::vector<NodeId> const& leaves() const {
        return _leaves;
    }

private:
    std::vector<std::vector<NodeId>> _adjacency_array;
    EdgeId                           _num_edges = 0;
    std::vector<NodeId>              _roots;
    std::vector<NodeId>              _leaves;
    std::vector<NodeId>              _preorder;
    std::vector<NodeId>              _postorder;
    std::vector<NodeId>              _levelorder;
};
