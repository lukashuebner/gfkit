#pragma once

// #include <sparsehash/dense_hash_map>
#include <memory>
#include <unordered_set>

#include <kassert/kassert.hpp>
#include <sfkit/include-redirects/cereal.hpp>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/checking_casts.hpp"
#include "sfkit/graph/AdjacencyArrayGraph.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/load/SubtreeHasher.hpp"
#include "sfkit/samples/NumSamplesBelow.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/tskit.hpp"
#include "sfkit/utils/concepts.hpp"

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

    [[nodiscard]] bool is_sample(NodeId const node_id) const {
        return _dag_postorder_edges.is_leaf(node_id);
    }

    [[nodiscard]] SampleSet all_samples() const {
        KASSERT(
            num_nodes_is_set(),
            "The number of nodes in the DAG must be set before calling all_samples().",
            sfkit::assert::light
        );
        SampleSet sample_set(num_samples());
        for (NodeId const sample_id: _dag_postorder_edges.leaves()) {
            KASSERT(sample_id < num_samples(), "Sample ID is out of range.", sfkit::assert::light);
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

    // TODO We're adding our leaves first, so they have ids [0,num_leaves). Use this to speed up the is_sample() query
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
    EdgeListGraph _dag_postorder_edges;
    // TODO Remove these (will cause a version bump in the  serialization format)
    std::vector<NodeId> _subtree_sizes;
};
