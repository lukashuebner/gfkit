#pragma once

#include <algorithm>
#include <cstddef>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <kassert/kassert.hpp>

#include "common.hpp"
#include "tdt/assertion_levels.hpp"

using EdgeList = std::vector<Edge>;

// TODO write unit tests
class EdgeListGraph {
public:
    using iterator       = EdgeList::iterator;
    using const_iterator = EdgeList::const_iterator;

    EdgeListGraph() {}

    EdgeListGraph(TraversalOrder traversal_order) : _traversal_order(traversal_order) {}

    void add_edge(NodeId from, NodeId to) {
        _edges.emplace_back(Edge(from, to));
        _num_nodes.reset();
    }

    void add_root(NodeId root) {
        _roots.push_back(root);
        _num_nodes.reset();
    }

    void add_leaf(NodeId leaf) {
        _leaves.push_back(leaf);
        _num_nodes.reset();
    }

    std::size_t num_edges() const {
        return _edges.size();
    }

    // TODO Cache the result of this computation.
    // TODO Write unit tests
    // Return the vertices in the graph as well as their outdegree.
    std::unordered_map<NodeId, std::size_t> nodes() const {
        std::unordered_map<NodeId, std::size_t> nodes;

        auto insert_or_increment = [&nodes](NodeId const vertex, bool increment = true) {
            auto&& it = nodes.find(vertex);
            if (it == nodes.end()) {
                nodes.insert({vertex, increment});
            } else if (increment) {
                it++;
            }
        };

        for (auto const& edge: _edges) {
            insert_or_increment(edge.from(), true);
            insert_or_increment(edge.to(), false);
        }

        for (auto const& root: _roots) {
            insert_or_increment(root, false);
        }

        for (auto const& leaf: _leaves) {
            insert_or_increment(leaf, false);
        }

        _num_nodes.emplace(nodes.size());
        return nodes;
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

    std::size_t num_roots() const {
        KASSERT(_unique_nodes(_roots), "Roots are not unique", tdt::assert::heavy);
        return _roots.size();
    }

    std::size_t num_trees() const {
        return num_roots();
    }

    std::vector<NodeId> const& leaves() const {
        KASSERT(_unique_nodes(_leaves), "Leaves are not unique", tdt::assert::heavy);
        return _leaves;
    }

    std::size_t num_leaves() const {
        KASSERT(_unique_nodes(_leaves), "Leaves are not unique", tdt::assert::heavy);
        return _leaves.size();
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

    void num_nodes(size_t num_nodes) {
        _num_nodes = num_nodes;
    }

    size_t num_nodes() const {
        if (_num_nodes.has_value()) {
            // The value is cached -> O(1)
            return _num_nodes.value();
        } else {
            // We have to recompute the node count -> O(|_edges| + |_leaves| + |_roots|)
            return nodes().size();
        }
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
    void serialize(Archive& ar) {
        ar(_edges, _roots, _leaves, _traversal_order);
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
    mutable std::optional<size_t> _num_nodes = std::nullopt; // This get's updated when nodes() is called.
    // TODO Cache the list of nodes?
    EdgeList            _edges;
    std::vector<NodeId> _roots;
    std::vector<NodeId> _leaves;
    TraversalOrder      _traversal_order = TraversalOrder::Unordered;
};
