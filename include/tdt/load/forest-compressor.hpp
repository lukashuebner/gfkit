#pragma once

#include <sparsehash/dense_hash_map>
#include <unordered_set>

#include <kassert/kassert.hpp>
// TODO Remove pasta::bit_set from the repo again

#include "tdt/assertion_levels.hpp"
#include "tdt/checking_casts.hpp"
#include "tdt/graph/adjacency-array-graph.hpp"
#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/load/subtree-id.hpp"
#include "tdt/tskit.hpp"
#include "tdt/utils/concepts.hpp"
#include "tdt/graph/compressed-forest.hpp"

// TODO Separate this class into the construction and the storage
class ForestCompressor {
    class SubtreeIdNodeMapper {
    public:
        SubtreeIdNodeMapper() {
            // TODO Try absl::flat_hash_map or https://github.com/Tessil/robin-map
            _subtree_to_node_map.set_empty_key(SuccinctSubtreeIdZero);
            // Even if I know the exact size of the map, reserving the memory /degrades/ performance.
            // Hypothesis: Even more cache-misses in the beginning, when the map isn't fully filled yet.
            // _subtree_to_node_map.resize(...);
        }

        NodeId insert(SuccinctSubtreeId const& subtree_id) {
            KASSERT(!contains(subtree_id), "Subtree ID already exists in the map", tdt::assert::light);
            _subtree_to_node_map[subtree_id] = _next_node_id;
            return _next_node_id++;
        }

        // Root nodes might have identical subtrees, thus the mapping is lo longer surjective. As roots are never
        // referred to, we do not need to store them but only assign them a node id.
        NodeId insert_root() {
            return _next_node_id++;
        }

        bool contains(SuccinctSubtreeId const& subtree_id) const {
            return _subtree_to_node_map.find(subtree_id) != _subtree_to_node_map.end();
        }

        NodeId get(SuccinctSubtreeId const& subtree_id) const {
            KASSERT(contains(subtree_id), "Subtree ID does not exists in the map", tdt::assert::light);
            return _subtree_to_node_map[subtree_id];
        }

        NodeId operator[](SuccinctSubtreeId const& subtree_id) const {
            return get(subtree_id);
        }

        NodeId num_nodes() const {
            return _next_node_id;
        }

    private:
        mutable google::dense_hash_map<SuccinctSubtreeId, NodeId> _subtree_to_node_map;
        NodeId                                                    _next_node_id = 0;
    };

public:
    ForestCompressor(TSKitTreeSequence& tree_sequence) : _tree_sequence(tree_sequence) {}

    CompressedForest compress() {
        CompressedForest         compressed_forest;
        TSKitTree                ts_tree(_tree_sequence);
        SuccinctSubtreeIdFactory succinct_subtree_id_factory;
        // Reserve the memory for the subtree id once. This prevents reallocations and re-filling the memory for every
        // subtree id, thus improves performance.
        // TODO Move this down for better readability, this optimization is not needed anymore.
        // Cache the (DAG) subtree IDs for this tree. As in tskit the same node can have different subtrees below it
        // depending on the current tree, we need to reset this map for each tree. (The subtree ids are in the DAG
        // domain in which each node id occurs only once)
        // TODO Write assertions for that no old subtree id is reused
        std::vector<SuccinctSubtreeId> ts_node_to_dag_subtree_id_map(ts_tree.max_node_id() + 1);

        _ts_node2cf_subtree.reserve(_tree_sequence.num_trees());
        bool first_tree = true;
        // TODO Rewrite this, once we have the tree_sequence iterator
        for (ts_tree.first(); ts_tree.is_valid(); ts_tree.next(), first_tree = false) {
            _ts_node2cf_subtree.resize(_ts_node2cf_subtree.size() + 1);
            _ts_node2cf_subtree.back().set_empty_key(TSK_NULL);
            // TODO Why is KASSERT not turned off completely in release mode?

            // TODO Decide on vertex vs node and use it consistently
            for (auto const ts_node_id: ts_tree.postorder()) {
                // Prefetching the children of the next node slows down the code in benchmarks.
                // TODO Use custom data structure for this, as the tskit one produces a lot of cache-misses
                if (ts_tree.is_sample(ts_node_id)) {
                    // Compute the subtree id of this leaf node by hashing its label.
                    // TODO Cache the result of this evaluation? There are a few thousand subtree IDs.
                    auto dag_subtree_id = succinct_subtree_id_factory.compute(ts_node_id);

                    // Cache the subtree id for this ts node
                    KASSERT(asserting_cast<size_t>(ts_node_id) < ts_node_to_dag_subtree_id_map.size());
                    ts_node_to_dag_subtree_id_map[asserting_cast<size_t>(ts_node_id)] = dag_subtree_id;

                    // Map the DAG subtree ID to the corresponding node ID in the DAG. The first tree should insert all
                    // the samples as the samples in all trees are identical.
                    if (first_tree) {
                        KASSERT(
                            !_dag_subtree_to_node_map.contains(dag_subtree_id),
                            "A leaf is present twice in the first tree.",
                            tdt::assert::normal
                        );
                        NodeId dag_node_id = _dag_subtree_to_node_map.insert(dag_subtree_id);
                        compressed_forest.add_leaf(dag_node_id);
                        _ts_node2cf_subtree.back()[ts_node_id] = dag_node_id;
                    } else {
                        KASSERT(
                            _dag_subtree_to_node_map.contains(dag_subtree_id),
                            "A leaf present in a later tree is missing from the first tree.",
                            tdt::assert::light
                        );
                        NodeId dag_node_id                     = _dag_subtree_to_node_map.get(dag_subtree_id);
                        _ts_node2cf_subtree.back()[ts_node_id] = dag_node_id;
                    }
                } else { // Node is inner node
                    // Compute the subtree id of this inner node by hashing the subtree ids of it's children.
                    // TODO is there a more elegant solution?
                    static_assert(std::is_same_v<decltype(XXH128_hash_t::low64), decltype(XXH128_hash_t::high64)>);
                    // std::vector<decltype(XXH128_hash_t::low64)> children_dag_subtree_ids;
                    SuccinctSubtreeId children_dag_subtree_ids = SuccinctSubtreeIdZero;
                    // TODO Don't reserve memory for this. Save the current position in the edge list, add the target
                    // nodes without the from edge and fix the from edge at the end of this block.
                    std::vector<NodeId> children_dag_node_ids;

                    // TODO Add checks for compression efficiency
                    // TODO Add compression of example input files
                    // TODO Prefetch children
                    for (auto child_ts_id: ts_tree.children(ts_node_id)) {
                        auto const child_dag_subtree_id =
                            ts_node_to_dag_subtree_id_map[asserting_cast<size_t>(child_ts_id)];

                        children_dag_subtree_ids ^= child_dag_subtree_id;

                        KASSERT(_dag_subtree_to_node_map.contains(child_dag_subtree_id));
                        children_dag_node_ids.emplace_back(_dag_subtree_to_node_map[child_dag_subtree_id]);
                    }
                    // TODO Instead of a full hash, maybe use only the hash's finisher. The values are random anyhow, as
                    // they were generated by xor-ing the hashes of the children.
                    auto dag_subtree_id = succinct_subtree_id_factory.compute(children_dag_subtree_ids);

                    // Map the DAG subtree ID to the node ID in the DAG.
                    KASSERT(asserting_cast<size_t>(ts_node_id) < ts_node_to_dag_subtree_id_map.size());
                    ts_node_to_dag_subtree_id_map[asserting_cast<size_t>(ts_node_id)] = dag_subtree_id;

                    // Add this node to the DAG if not already present. As the DAG is stored
                    // as a list of edges, we need to add an edge from this node to each of
                    // it's children.  In the case that two trees in the tree sequence are exactly identical, we want
                    // two root nodes in the DAG -- one for each of the two trees.
                    bool const is_root = ts_tree.is_root(ts_node_id);
                    // TODO We are querying the hash map twice here. Rewrite code to not do this.
                    if (!_dag_subtree_to_node_map.contains(dag_subtree_id) || is_root) {
                        NodeId dag_node_id;
                        // Root nodes can have the same ID (if both trees are identical), but don't have in edges. Thus,
                        // we don't need to map their ID.
                        if (is_root) {
                            dag_node_id = _dag_subtree_to_node_map.insert_root();
                            // If the node is a root node in the tree sequence, also add it to the DAG as a root node.
                            compressed_forest.add_root(dag_node_id);
                        } else {
                            dag_node_id = _dag_subtree_to_node_map.insert(dag_subtree_id);
                        }
                        for (auto child_dag_node_id: children_dag_node_ids) {
                            compressed_forest.add_edge(dag_node_id, child_dag_node_id);
                        }
                        _ts_node2cf_subtree.back()[ts_node_id] = dag_node_id;
                    } else {
                        _ts_node2cf_subtree.back()[ts_node_id] = _dag_subtree_to_node_map[dag_subtree_id];
                    }
                }
            }
        }
        KASSERT(compressed_forest.roots().size() == _tree_sequence.num_trees());
        KASSERT(compressed_forest.num_roots() == _tree_sequence.num_trees());
        KASSERT(compressed_forest.num_leaves() == _tree_sequence.num_samples());

        // Set the number of nodes in the DAG so it does not have to be recomputed.
        compressed_forest.num_nodes(_dag_subtree_to_node_map.num_nodes());

        // As we build the tree edges by a postorder traversal on the tree, the from edges should be post-ordered, too.
        compressed_forest.postorder_edges().traversal_order(TraversalOrder::Postorder);

        return compressed_forest;
    }

    // TODO Rework the naming of nodes and subtrees in the tree sequence and the dag
    NodeId ts_node2cf_subtree(TreeId tree_id, tsk_id_t ts_node) const {
        return _ts_node2cf_subtree[asserting_cast<size_t>(tree_id)][ts_node];
    }

    template <IterableInput IterableInput>
    std::vector<NodeId> ts_node2cf_subtree(TreeId tree_id, IterableInput input) const {
        std::vector<NodeId> result;
        result.reserve(input.size());
        std::transform(input.begin(), input.end(), std::back_inserter(result), [this, tree_id](auto&& tsk_node) {
            return ts_node2cf_subtree(tree_id, tsk_node);
        });
        return result;
    }

private:
    SubtreeIdNodeMapper _dag_subtree_to_node_map;
    TSKitTreeSequence&  _tree_sequence;

    // TODO Move this to a separate class
    // TODO Use a different data structure for this?
    mutable std::vector<google::dense_hash_map<tsk_id_t, NodeId>> _ts_node2cf_subtree;
};
