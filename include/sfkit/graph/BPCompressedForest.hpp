#pragma once

// #include <sparsehash/dense_hash_map>
#include <memory>
#include <unordered_set>

#include <kassert/kassert.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sfkit/include-redirects/cereal.hpp>

// TODO include-what-you-use
#include "sfkit/assertion_levels.hpp"
#include "sfkit/checking_casts.hpp"
#include "sfkit/graph/AdjacencyArrayGraph.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/load/SubtreeHashToNodeMapper.hpp"
#include "sfkit/load/SubtreeHasher.hpp"
#include "sfkit/samples/NumSamplesBelow.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/tskit.hpp"
#include "sfkit/utils/concepts.hpp"

class BPCompressedForest {
public:
    // TODO We don't move these vectors as we want to compress them. -> Add compression
    // TODO add const where possible
    BPCompressedForest(
        sdsl::bit_vector const&       is_reference,
        sdsl::bit_vector const&       balanced_parenthesis,
        sdsl::int_vector<> const&     references,
        sdsl::int_vector<> const&     leaves,
        SubtreeHashToNodeMapper const& subtree_hash_to_node_mapper,
        NodeId const                  num_nodes,
        NodeId const                  num_leaves
    )
        : _is_reference(is_reference),
          _balanced_parenthesis(balanced_parenthesis),
          _references(references),
          _leaves(leaves),
          _num_nodes(num_nodes),
          _num_leaves(num_leaves) {}

    [[nodiscard]] NodeId num_nodes() const {
        return _num_nodes;
    }

    // void num_nodes(NodeId num_nodes) {
    //     return _dag_postorder_edges.num_nodes(num_nodes);
    // }

    // [[nodiscard]] bool num_nodes_is_set() const {
    //     return _dag_postorder_edges.num_nodes_is_set();
    // }

    // [[nodiscard]] bool is_sample(NodeId const node_id) const {
    //     return _dag_postorder_edges.is_leaf(node_id);
    // }

    [[nodiscard]] SampleSet all_samples() const {
        SampleSet sample_set(num_samples());
        // TODO Re-add
        // for (NodeId const sample_id: _leaves) {
        //     KASSERT(sample_id < num_samples(), "Sample ID is out of range.", sfkit::assert::light);
        //     sample_set.add(sample_id);
        // }
        return sample_set;
    }

    [[nodiscard]] SampleId num_samples() const {
        return _num_leaves;
    }

    // [[nodiscard]] TreeId num_trees() const {
    //     return _dag_postorder_edges.num_trees();
    // }

    // [[nodiscard]] EdgeId num_edges() const {
    //     return _dag_postorder_edges.num_edges();
    // }

    // // TODO We're adding our leaves first, so they have ids [0,num_leaves). Use this to speed up the is_sample()
    // query void insert_leaf(NodeId const leaf) {
    //     _dag_postorder_edges.insert_leaf(leaf);
    // }

    // void insert_root(NodeId const root) {
    //     _dag_postorder_edges.insert_root(root);
    // }

    // void insert_edge(NodeId const from, NodeId const to) {
    //     _dag_postorder_edges.insert_edge(from, to);
    // }

    // [[nodiscard]] std::vector<NodeId> const& roots() const {
    //     return _dag_postorder_edges.roots();
    // }

    // [[nodiscard]] NodeId num_roots() const {
    //     return _dag_postorder_edges.num_roots();
    // }

    // [[nodiscard]] std::vector<NodeId> const& leaves() const {
    //     return _dag_postorder_edges.leaves();
    // }

    [[nodiscard]] SampleId num_leaves() const {
        return _num_leaves;
    }

    [[nodiscard]] NodeId num_unique_subtrees() const {
        return _num_nodes;
    }

    // [[nodiscard]] std::unordered_set<NodeId> unique_subtrees() const {
    //     return _dag_postorder_edges.nodes();
    // }

    // This function is accessible mainly for unit-testing. It is not part of the public API.
    auto const is_reference() const {
        return _is_reference;
    }

    auto const balanced_parenthesis() const {
        return _balanced_parenthesis;
    }

    auto const references() const {
        return _references;
    }

    auto const leaves() const {
        return _leaves;
    }

private:
    sdsl::bit_vector   _is_reference;
    sdsl::bit_vector   _balanced_parenthesis;
    sdsl::int_vector<> _references;
    sdsl::int_vector<> _leaves;
    size_t             _num_nodes;
    size_t             _num_leaves;
};
