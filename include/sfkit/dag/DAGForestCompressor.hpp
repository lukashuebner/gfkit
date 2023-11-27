#pragma once

// #include <sparsehash/dense_hash_map>
#include <unordered_set>

#include <kassert/kassert.hpp>
#include <sfkit/include-redirects/hopscotch_map.hpp>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/dag/DAGCompressedForest.hpp"
#include "sfkit/graph/AdjacencyArrayGraph.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/graph/SubtreeHashToNodeMapper.hpp"
#include "sfkit/graph/SubtreeHasher.hpp"
#include "sfkit/graph/TsToSfNodeMapper.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"
#include "sfkit/sequence/GenomicSequenceFactory.hpp"
#include "sfkit/tskit/tskit.hpp"
#include "sfkit/utils/checking_casts.hpp"
#include "sfkit/utils/concepts.hpp"

namespace sfkit::dag {

using sfkit::utils::asserting_cast;
using namespace sfkit::graph;

class DAGForestCompressor {
public:
    DAGForestCompressor(tskit::TSKitTreeSequence& tree_sequence)
        : _tree_sequence(tree_sequence),
          _num_samples(tree_sequence.num_samples()),
          _ts_tree(tree_sequence) {
        if (!tree_sequence.sample_ids_are_consecutive()) {
            throw std::runtime_error("Sample IDs of the tree sequence are not consecutive.");
        }
        _ts_node_to_subtree.resize(_ts_tree.max_node_id());
    }

    template <typename GenomicSequenceFactoryT>
    DAGCompressedForest compress(GenomicSequenceFactoryT& genomic_sequence_factory) {
        DAGCompressedForest forest;

        // Add the samples to the compressed forest here so they have the same ID in all trees.
        _register_samples(forest);

        // TODO Rewrite this, once we have the tree_sequence iterator
        for (_ts_tree.first(); _ts_tree.is_valid(); _ts_tree.next()) {
            for (auto const ts_node_id: _ts_tree.postorder()) {
                // Samples are already mapped and added to the DAG before processing the first tree.
                if (is_sample(ts_node_id)) [[unlikely]] {
                    continue;
                }

                // Compute the subtree ID of this inner node by hashing the subtree IDs of its children.
                _subtree_hash_factory.reset();

                // TODO Don't reserve memory for this. Save the current position in the edge list, add the target
                // nodes without the from edge and fix the from edge at the end of this block.
                std::vector<NodeId> children_sf_node_ids;

                for (auto child_ts_id: _ts_tree.children(ts_node_id)) {
                    auto const childs_subtree_id = _ts_node_to_subtree[asserting_cast<size_t>(child_ts_id)];

                    _subtree_hash_factory.append_child(childs_subtree_id);

                    KASSERT(_subtree_to_sf_node.contains(childs_subtree_id));
                    children_sf_node_ids.emplace_back(_subtree_to_sf_node[childs_subtree_id]);
                }
                SubtreeHash const subtree_id = _subtree_hash_factory.hash();

                // Map the DAG subtree ID to the node ID in the DAG.
                KASSERT(asserting_cast<size_t>(ts_node_id) < _ts_node_to_subtree.size());
                _ts_node_to_subtree[asserting_cast<size_t>(ts_node_id)] = subtree_id;

                // Add this node to the DAG if not already present. As the DAG is stored
                // as a list of edges, we need to add an edge from this node to each of
                // its children.  In the case that two trees in the tree sequence are exactly identical, we want
                // two root nodes in the DAG -- one for each of the two trees.
                bool const subtree_is_root = _ts_tree.is_root(ts_node_id);
                auto const sf_node_it      = _subtree_to_sf_node.find(subtree_id);
                bool const subtree_in_dag  = sf_node_it != _subtree_to_sf_node.end();

                // If the subtree is not in the DAG or if it is a root node, add it to the DAG.
                if (!subtree_in_dag || subtree_is_root) [[unlikely]] {
                    NodeId sf_node_id = INVALID_NODE_ID;

                    // Root nodes can have the same ID (if both trees are identical), but don't have in edges. Thus,
                    // we don't need to map their ID.
                    if (subtree_is_root) [[unlikely]] {
                        sf_node_id = _subtree_to_sf_node.insert_or_update_node(subtree_id);

                        // If the node is a root node in the tree sequence, also add it to the DAG as a root node.
                        forest.insert_root(sf_node_id);
                    } else {
                        sf_node_id = _subtree_to_sf_node.insert_node(subtree_id);
                    }
                    KASSERT(sf_node_id != INVALID_NODE_ID);

                    // Add the edges from this sf node to its children's sf nodes to the DAG
                    for (auto&& child: children_sf_node_ids) {
                        forest.insert_edge(sf_node_id, child);
                    }
                }
            }

            // Process the mutations of this tree
            genomic_sequence_factory.process_mutations(
                asserting_cast<TreeId>(_ts_tree.tree_id()),
                graph::TsToSfNodeMapper(_ts_node_to_subtree, _subtree_to_sf_node)
            );
        }

        genomic_sequence_factory.finalize();

        KASSERT(forest.roots().size() == _tree_sequence.num_trees());
        KASSERT(forest.num_roots() == _tree_sequence.num_trees());
        KASSERT(forest.num_leaves() == _tree_sequence.num_samples());

        // Set the number of nodes in the DAG so it does not have to be recomputed.
        forest.num_nodes(_subtree_to_sf_node.num_nodes());

        // As we build the tree edges by a postorder traversal on the tree, the from edges should be post-ordered, too.
        forest.postorder_edges().traversal_order(TraversalOrder::Postorder);

        return forest;
    }

private:
    tskit::TSKitTreeSequence& _tree_sequence;
    tsk_size_t                _num_samples;
    tskit::TSKitTree          _ts_tree;
    std::vector<SubtreeHash>  _ts_node_to_subtree;
    SubtreeHashToNodeMapper   _subtree_to_sf_node;
    SubtreeHasher             _subtree_hash_factory;

    // The sample ids are consecutive: 0 ... num_samples - 1
    inline bool is_sample(tsk_id_t ts_node_id) const {
        return ts_node_id < asserting_cast<tsk_id_t>(_num_samples);
    }

    // Add them to the compressed forest first, so they have the same IDs there.
    void _register_samples(DAGCompressedForest& forest) {
        KASSERT(_tree_sequence.sample_ids_are_consecutive(), "Sample IDs are not consecutive.");
        for (SampleId sample_id = 0; sample_id < _num_samples; sample_id++) {
            KASSERT(_tree_sequence.is_sample(asserting_cast<tsk_id_t>(sample_id)));

            // Compute the subtree ID of this sample (leaf) node by hashing its label.
            auto subtree_hash = _subtree_hash_factory.hash_sample(sample_id);

            // Cache the subtree ID for this TS node
            KASSERT(asserting_cast<size_t>(sample_id) < _ts_node_to_subtree.size());
            _ts_node_to_subtree[asserting_cast<size_t>(sample_id)] = subtree_hash;

            // Map the DAG subtree ID to the corresponding node ID in the DAG. The first tree should insert all
            // the samples as the samples in all trees are identical.
            KASSERT(
                !_subtree_to_sf_node.contains(subtree_hash),
                "A leaf is present twice in the first tree.",
                sfkit::assert::normal
            );
            NodeId dag_node_id = _subtree_to_sf_node.insert_node(subtree_hash);
            forest.insert_leaf(dag_node_id);
        }
    }
};
} // namespace sfkit::dag
