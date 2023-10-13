#pragma once

// #include <sparsehash/dense_hash_map>
#include <memory>
#include <unordered_set>

#if defined(__GNUC__) && !defined(__clang__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wnoexcept"
#endif
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#if defined(__GNUC__) && !defined(__clang__)
    #pragma GCC diagnostic pop
#endif
#include <kassert/kassert.hpp>

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

    void num_nodes(NodeId num_nodes) {
        return _dag_postorder_edges.num_nodes(num_nodes);
    }

    [[nodiscard]] bool num_nodes_is_set() const {
        return _dag_postorder_edges.num_nodes_is_set();
    }

    // [[nodiscard]] std::shared_ptr<NumSamplesBelow> compute_num_samples_below(SampleSet const& sample_set) const {
    //     return std::make_shared<NumSamplesBelow>(_dag_postorder_edges, sample_set);
    // }

    [[nodiscard]] bool is_sample(NodeId const node_id) const {
        return _dag_postorder_edges.is_leaf(node_id);
    }

    [[nodiscard]] SampleSet all_samples() const {
        KASSERT(
            num_nodes_is_set(),
            "The number of nodes in the DAG must be set before calling all_samples().",
            tdt::assert::light
        );
        SampleSet sample_set(num_samples());
        for (NodeId const sample_id: _dag_postorder_edges.leaves()) {
            KASSERT(sample_id < num_samples(), "Sample ID is out of range.", tdt::assert::light);
            sample_set.add(sample_id);
        }
        return sample_set;
    }

    [[nodiscard]] SampleId num_samples() const {
        return _dag_postorder_edges.num_leaves();
    }

    [[nodiscard]] TreeId num_trees() const {
        return _dag_postorder_edges.num_trees();
    }

    [[nodiscard]] EdgeId num_edges() const {
        return _dag_postorder_edges.num_edges();
    }

    void insert_leaf(NodeId const leaf) {
        _dag_postorder_edges.insert_leaf(leaf);
    }

    void insert_root(NodeId const root) {
        _dag_postorder_edges.insert_root(root);
    }

    void insert_edge(NodeId const from, NodeId const to) {
        _dag_postorder_edges.insert_edge(from, to);
    }

    [[nodiscard]] std::vector<NodeId> const& roots() const {
        return _dag_postorder_edges.roots();
    }

    [[nodiscard]] NodeId num_roots() const {
        return _dag_postorder_edges.num_roots();
    }

    [[nodiscard]] std::vector<NodeId> const& leaves() const {
        return _dag_postorder_edges.leaves();
    }

    [[nodiscard]] SampleId num_leaves() const {
        return _dag_postorder_edges.num_leaves();
    }

    [[nodiscard]] NodeId num_unique_subtrees() const {
        return _dag_postorder_edges.num_nodes();
    }

    [[nodiscard]] std::unordered_set<NodeId> unique_subtrees() const {
        return _dag_postorder_edges.nodes();
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
