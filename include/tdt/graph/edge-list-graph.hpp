#pragma once

#include <algorithm>
#include <cstddef>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include <tsl/hopscotch_map.h>
#pragma GCC diagnostic pop

#include <kassert/kassert.hpp>

#include "common.hpp"
#include "tdt/assertion_levels.hpp"
#include "tdt/samples/sample-set.hpp"

using EdgeList = std::vector<Edge>;

// TODO write unit tests
class EdgeListGraph {
public:
    using iterator         = EdgeList::iterator;
    using const_iterator   = EdgeList::const_iterator;
    using NodeOutdegreeMap = tsl::hopscotch_map<NodeId, NodeId>;

    EdgeListGraph() {}

    EdgeListGraph(TraversalOrder traversal_order) : _traversal_order(traversal_order) {}

    void insert_edge(NodeId from, NodeId to) {
        _edges.emplace_back(Edge(from, to));
        // TODO Think about incremental updates of the nodes array or a separate build phase for the graph.
    }

    void insert_root(NodeId root) {
        _roots.push_back(root);
    }

    void insert_leaf(NodeId leaf) {
        _leaves.push_back(leaf);
    }

    EdgeId num_edges() const {
        return asserting_cast<EdgeId>(_edges.size());
    }

    iterator begin() noexcept {
        return _edges.begin();
    }

    iterator end() noexcept {
        return _edges.end();
    }

    const_iterator begin() const noexcept {
        return _edges.begin();
    }

    const_iterator end() const noexcept {
        return _edges.end();
    }

    const_iterator cbegin() const noexcept {
        return _edges.cbegin();
    }

    const_iterator cend() const noexcept {
        return _edges.cend();
    }

    std::vector<NodeId> const& roots() const {
        KASSERT(_unique_nodes(_roots), "Roots are not unique", tdt::assert::heavy);
        return _roots;
    }

    NodeId num_roots() const {
        KASSERT(_unique_nodes(_roots), "Roots are not unique", tdt::assert::heavy);
        return asserting_cast<NodeId>(_roots.size());
    }

    TreeId num_trees() const {
        return num_roots();
    }

    std::vector<NodeId> const& leaves() const {
        KASSERT(_unique_nodes(_leaves), "Leaves are not unique", tdt::assert::heavy);
        return _leaves;
    }

    SampleId num_leaves() const {
        KASSERT(_unique_nodes(_leaves), "Leaves are not unique", tdt::assert::heavy);
        return asserting_cast<SampleId>(_leaves.size());
    }

    bool directed() const {
        return true;
    }

    TraversalOrder traversal_order() const {
        return _traversal_order;
    }

    void traversal_order(TraversalOrder traversal_order) {
        _traversal_order = traversal_order;
    }

    bool is_inorder() const {
        return _traversal_order == TraversalOrder::Inorder;
    }

    bool is_postorder() const {
        return _traversal_order == TraversalOrder::Postorder;
    }

    bool is_preorder() const {
        return _traversal_order == TraversalOrder::Preorder;
    }

    bool is_levelorder() const {
        return _traversal_order == TraversalOrder::Levelorder;
    }

    bool is_unordered() const {
        return _traversal_order == TraversalOrder::Unordered;
    }

    enum class SortBy {
        FromVertex,
        ToVertex,
    };

    // TODO make this const
    std::unordered_set<NodeId> nodes() const {
        std::unordered_set<NodeId> nodes;
        nodes.reserve(2 * _leaves.size() + (_roots.size() - 1)); // Lower bound only

        for (auto const& edge: _edges) {
            nodes.insert(edge.from());
            nodes.insert(edge.to());
        }

        for (auto const& root: _roots) {
            nodes.insert(root);
        }

        for (auto const& leaf: _leaves) {
            nodes.insert(leaf);
        }

        return nodes;
    }

    // We have to recompute the node count -> O(|_edges| + |_leaves| + |_roots|)
    // This is incredibly slow, so use it only in unit tests
    void compute_num_nodes() {
        _num_nodes = asserting_cast<NodeId>(nodes().size());
    }

    void num_nodes(NodeId num_nodes) {
        KASSERT(!num_nodes_is_set(), "The number of nodes is already set", tdt::assert::light);
        _num_nodes = num_nodes;
    }

    [[nodiscard]] NodeId num_nodes() const {
        KASSERT(num_nodes_is_set(), "The number of nodes is not set", tdt::assert::light);
        return _num_nodes;
    }

    [[nodiscard]] bool num_nodes_is_set() const {
        return _num_nodes != INVALID_NODE_ID;
    }

    bool check_postorder() const {
        // Initialize all leafs as visited and all other nodes as unvisited.
        std::vector<bool> visited(num_nodes(), false);
        for (auto leaf: _leaves) {
            visited[leaf] = true;
        }

        // Traverse the graph in the order given by _edges
        for (auto const& edge: _edges) {
            // In postorder, the parent node must be visited /after/ the child node.
            if (!visited[edge.to()]) {
                return false;
            }
            visited[edge.from()] = true;
        }

        // Check that all roots were visited.
        for (auto root: _roots) {
            if (!visited[root]) {
                return false;
            } else {
                visited[root] = true; // Mark the root as visited for the check below.
            }
        }

        // Check that all nodes were visited -- this does not check if the leafs are visited, as we set them to true
        // explicitly above.
        for (auto const& v: visited) {
            if (!v) {
                return false;
            }
        }

        return true;
    }

    template <class Compare>
    void sort_edges(Compare comp, TraversalOrder traversal_order = TraversalOrder::Unordered) {
        _traversal_order = traversal_order;
        std::sort(_edges.begin(), _edges.end(), comp);
    }

    void sort_edges(SortBy sort_by) {
        switch (sort_by) {
            case SortBy::FromVertex:
                sort_edges([](Edge const& a, Edge const& b) { return a.from() < b.from(); }, TraversalOrder::Unordered);
                break;
            case SortBy::ToVertex:
                sort_edges([](Edge const& a, Edge const& b) { return a.to() < b.to(); }, TraversalOrder::Unordered);
                break;
            default:
                KASSERT(false, "Invalid sort_by value", tdt::assert::light);
                break;
        }
    }

    bool edges_are_sorted(SortBy sort_by) const {
        switch (sort_by) {
            case SortBy::FromVertex:
                return std::is_sorted(_edges.begin(), _edges.end(), [](Edge const& a, Edge const& b) {
                    return a.from() < b.from();
                });
                break;
            case SortBy::ToVertex:
                return std::is_sorted(_edges.begin(), _edges.end(), [](Edge const& a, Edge const& b) {
                    return a.to() < b.to();
                });
                break;
            default:
                KASSERT(false, "Invalid sort_by value", tdt::assert::light);
                return false;
                break;
        }
    }

    bool is_leaf(NodeId const node) const {
        return std::find(_leaves.begin(), _leaves.end(), node) != _leaves.end();
    }

    template <class Archive>
    void serialize(Archive& archive) {
        archive(_num_nodes, _edges, _roots, _leaves, _traversal_order);
    }

private:
    // We don't use an inplace approach in order not to have a side effect only if assertions are enabled.
    bool _unique_nodes(std::vector<NodeId> const& nodes) const {
        std::unordered_set<NodeId> set;
        for (auto const& node: nodes) {
            if (set.contains(node)) {
                return false;
            }
            set.insert(node);
        }
        KASSERT(set.size() == nodes.size());
        return true;
    }

    // We're not assuming that node ids are consecutive
    NodeId              _num_nodes = INVALID_NODE_ID;
    EdgeList            _edges;
    std::vector<NodeId> _roots;
    std::vector<NodeId> _leaves;
    TraversalOrder      _traversal_order = TraversalOrder::Unordered;
};
