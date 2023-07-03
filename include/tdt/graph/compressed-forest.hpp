#pragma once

// #include <sparsehash/dense_hash_map>
#include <unordered_set>

#include <kassert/kassert.hpp>

#include "edge-list-graph.hpp"
#include "tdt/assertion_levels.hpp"
#include "tdt/checking_casts.hpp"
#include "tdt/graph/adjacency-array-graph.hpp"
#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/load/subtree-id.hpp"
#include "tdt/tskit.hpp"
#include "tdt/utils/concepts.hpp"

#include <cereal/types/vector.hpp>
#include <cereal/types/memory.hpp>

class CompressedForest {
public:
    EdgeListGraph const& postorder_edges() const {
        return _dag_postorder_edges;
    }

    EdgeListGraph& postorder_edges() {
        return _dag_postorder_edges;
    }

    void compute_num_samples_below() {
        KASSERT(_dag_postorder_edges.check_postorder(), "DAG edges are not post-ordered.", tdt::assert::normal);
        std::fill(_subtree_sizes.begin(), _subtree_sizes.end(), 0);
        _subtree_sizes.resize(_dag_postorder_edges.num_nodes(), 0);
        for (auto&& leaf: _dag_postorder_edges.leaves()) {
            _subtree_sizes[leaf] = 1;
        }
        // TODO Try prefetching a few edges ahead
        for (auto&& edge: _dag_postorder_edges) {
            _subtree_sizes[edge.from()] += _subtree_sizes[edge.to()];
            KASSERT(
                _subtree_sizes[edge.from()] <= _dag_postorder_edges.num_leaves(),
                "Number of samples below a node exceeds the number of samples in the tree sequence.",
                tdt::assert::light
            );
        }
    }

    size_t num_nodes() const {
        return _dag_postorder_edges.num_nodes();
    }

    void num_nodes(size_t const num_nodes) {
        _dag_postorder_edges.num_nodes(num_nodes);
    }

    size_t num_samples_below(SubtreeId subtree_id) {
        if (_subtree_sizes.size() == 0) {
            this->compute_num_samples_below();
        }
        KASSERT(subtree_id < _subtree_sizes.size(), "Subtree ID out of bounds.", tdt::assert::light);
        return _subtree_sizes[subtree_id];
    }

    auto& num_samples_below() {
        if (_subtree_sizes.size() == 0) {
            this->compute_num_samples_below();
        }
        return _subtree_sizes;
    }

    bool is_sample(NodeId const node_id) const {
        return _dag_postorder_edges.is_leaf(node_id);
    }

    size_t num_samples() const {
        return _dag_postorder_edges.num_leaves();
    }

    size_t num_trees() const {
        return _dag_postorder_edges.num_trees();
    }

    size_t num_edges() const {
        return _dag_postorder_edges.num_edges();
    }

    void add_leaf(NodeId const leaf) {
        _dag_postorder_edges.add_leaf(leaf);
    }

    void add_root(NodeId const root) {
        _dag_postorder_edges.add_root(root);
    }

    void add_edge(NodeId const from, NodeId const to) {
        _dag_postorder_edges.add_edge(from, to);
    }

    std::vector<NodeId> const& roots() const {
        return _dag_postorder_edges.roots();
    }

    std::size_t num_roots() const {
        return _dag_postorder_edges.num_roots();
    }

    std::vector<NodeId> const& leaves() const {
        return _dag_postorder_edges.leaves();
    }

    std::size_t num_leaves() const {
        return _dag_postorder_edges.num_leaves();
    }

    template <class Archive>
    void serialize(Archive& ar) {
        ar(_dag_postorder_edges, _subtree_sizes);
    }

private:
    EdgeListGraph       _dag_postorder_edges;
    std::vector<size_t> _subtree_sizes;
};
