#pragma once

// #include <sparsehash/dense_hash_map>
#include <unordered_set>

#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <kassert/kassert.hpp>

#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/assertion_levels.hpp"
#include "tdt/checking_casts.hpp"
#include "tdt/graph/adjacency-array-graph.hpp"
#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/load/subtree-id.hpp"
#include "tdt/samples/num-samples-below.hpp"
#include "tdt/samples/sample-set.hpp"
#include "tdt/tskit.hpp"
#include "tdt/utils/concepts.hpp"

class CompressedForest {
public:
    EdgeListGraph const& postorder_edges() const {
        return _dag_postorder_edges;
    }

    EdgeListGraph& postorder_edges() {
        return _dag_postorder_edges;
    }

    [[nodiscard]] NodeId num_nodes() const {
        return _dag_postorder_edges.num_nodes();
    }

    [[nodiscard]] bool nodes_are_computed() const {
        return _dag_postorder_edges.nodes_are_computed();
    }

    void compute_nodes() const {
        _dag_postorder_edges.compute_nodes();
    }

    NumSamplesBelow compute_num_samples_below(SampleSet const& sample_set) const {
        return NumSamplesBelow(_dag_postorder_edges, sample_set);
    }

    bool is_sample(NodeId const node_id) const {
        return _dag_postorder_edges.is_leaf(node_id);
    }

    SampleSet all_samples() const {
        SampleSet sample_set(num_nodes());
        for (NodeId const node_id : _dag_postorder_edges.leaves()) {
            sample_set.add(node_id);
        }
        return sample_set;
    }

    SampleId num_samples() const {
        return _dag_postorder_edges.num_leaves();
    }

    TreeId num_trees() const {
        return _dag_postorder_edges.num_trees();
    }

    EdgeId num_edges() const {
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

    NodeId num_roots() const {
        return _dag_postorder_edges.num_roots();
    }

    std::vector<NodeId> const& leaves() const {
        return _dag_postorder_edges.leaves();
    }

    SampleId num_leaves() const {
        return _dag_postorder_edges.num_leaves();
    }

    template <class Archive>
    void serialize(Archive& ar) {
        // The number of nodes in the DAG are computed during serialization of the EdgeListGraph object.
        ar(_dag_postorder_edges, _subtree_sizes);
    }

private:
    EdgeListGraph       _dag_postorder_edges;
    std::vector<NodeId> _subtree_sizes;
};
