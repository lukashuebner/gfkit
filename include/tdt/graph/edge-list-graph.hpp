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

    void add_edge(NodeId from, NodeId to) {
        _edges.emplace_back(Edge(from, to));
        // TODO Think about incremental updates of the nodes array or a separate build phase for the graph.
        _nodes_are_valid = false;
    }

    void add_root(NodeId root) {
        _roots.push_back(root);
        _nodes_are_valid = false;
    }

    void add_leaf(NodeId leaf) {
        _leaves.push_back(leaf);
        _nodes_are_valid = false;
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

    void compute_nodes() const {
        _nodes.clear();

        auto insert_or_increment = [this](NodeId const vertex, bool increment = true) {
            auto&& it = _nodes.find(vertex);
            if (it == _nodes.end()) {
                _nodes.insert({vertex, increment});
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

        _nodes_are_valid = true;
    }

    // TODO Write unit tests
    // Return the vertices in the graph as well as their outdegree.
    NodeOutdegreeMap nodes() const {
        KASSERT(_nodes_are_valid, "The nodes are not computed yet", tdt::assert::light);
        return _nodes;
    }

    // TODO introduce finalize step which computes this
    NodeId num_nodes() const {
        // We have to recompute the node count -> O(|_edges| + |_leaves| + |_roots|)
        KASSERT(_nodes_are_valid, "The nodes are not computed yet", tdt::assert::light);
        return asserting_cast<NodeId>(_nodes.size());
    }

    // TODO The nodes might not be up to date. The proper way to solve this is to introduce a dirty flag.
    [[nodiscard]] bool nodes_are_computed() const {
        return !_nodes.empty();
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
    void save(Archive& archive) const {
        // Serialization of tsl::hopscotch_map is not trivial, thus we do not add it to the archive and recompute the
        // nodes when loading instead.
        archive(_edges, _roots, _leaves, _traversal_order);
    }

    template <class Archive>
    void load(Archive& archive) {
        archive(_edges, _roots, _leaves, _traversal_order);
        compute_nodes();
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

    // TODO Should we? Would this provide a speedup?
    // We're not assuming that node ids are consecutive
    // TODO remove this mutable
    mutable NodeOutdegreeMap _nodes;
    mutable bool             _nodes_are_valid = true;
    EdgeList                 _edges;
    std::vector<NodeId>      _roots;
    std::vector<NodeId>      _leaves;
    TraversalOrder           _traversal_order = TraversalOrder::Unordered;
};
