#pragma once

// #include <sparsehash/dense_hash_map>
#include <unordered_set>

#include <kassert/kassert.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include <tsl/hopscotch_map.h>
#pragma GCC diagnostic pop

#include "tdt/assertion_levels.hpp"
#include "tdt/checking_casts.hpp"
#include "tdt/graph/adjacency-array-graph.hpp"
#include "tdt/graph/compressed-forest.hpp"
#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/load/subtree-id-2-node-mapper.hpp"
#include "tdt/load/subtree-id.hpp"
#include "tdt/load/ts-node-2-sf-subtree-mapper.hpp"
#include "tdt/sequence/genomic-sequence-storage-factory.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"
#include "tdt/tskit.hpp"
#include "tdt/utils/concepts.hpp"

// TODO Separate this class into the construction and the storage
class ForestCompressor {
public:
    ForestCompressor(TSKitTreeSequence& tree_sequence) : _tree_sequence(tree_sequence) {
        if (!tree_sequence.sample_ids_are_consecutive()) {
            throw std::runtime_error("Sample IDs of the tree sequence are not consecutive.");
        }
    }

    // TODO Decide on vertex vs node and use it consistently
    CompressedForest compress(GenomicSequenceFactory& genomic_sequence_storage_factory) {
        CompressedForest         compressed_forest;
        TSKitTree                ts_tree(_tree_sequence);
        SuccinctSubtreeIdFactory succinct_subtree_id_factory;
        // Reserve the memory for the subtree id once. This prevents reallocations and re-filling the memory for every
        // subtree id, thus improves performance.
        // TODO Write assertions for that no old subtree id is reused
        std::vector<SuccinctSubtreeId> ts_node_to_dag_subtree_id_map(ts_tree.max_node_id() + 1);

        // Cache the (DAG) subtree IDs for the current tree. As in tskit the same node can have different subtrees below
        // it depending on the current tree, we need to reset this map for each tree. (The subtree ids are in the DAG
        // domain in which each node id occurs only once)
        TsNode2SfSubtreeMapper ts_node2sf_subtree;
        // TODO reserve space?

        std::vector<SubtreeId> samples_subtree_ids(_tree_sequence.num_samples());

        // The sample ids are consecutive: 0 ... num_samples - 1
        // Add them to the compressed forest first, so they have the same IDs there.
        KASSERT(_tree_sequence.sample_ids_are_consecutive(), "Sample IDs are not consecutive.");
        for (SampleId sample_id = 0; sample_id < _tree_sequence.num_samples(); sample_id++) {
            KASSERT(_tree_sequence.is_sample(asserting_cast<tsk_id_t>(sample_id)));

            // Compute the subtree id of this leaf node by hashing its label.
            auto dag_subtree_id = succinct_subtree_id_factory.compute(sample_id);

            // Cache the subtree id for this ts node
            KASSERT(asserting_cast<size_t>(sample_id) < ts_node_to_dag_subtree_id_map.size());
            ts_node_to_dag_subtree_id_map[asserting_cast<size_t>(sample_id)] = dag_subtree_id;

            // Map the DAG subtree ID to the corresponding node ID in the DAG. The first tree should insert all
            // the samples as the samples in all trees are identical.
            KASSERT(
                !_dag_subtree_to_node_map.contains(dag_subtree_id),
                "A leaf is present twice in the first tree.",
                tdt::assert::normal
            );
            NodeId dag_node_id = _dag_subtree_to_node_map.insert(dag_subtree_id);
            compressed_forest.add_leaf(dag_node_id);
            ts_node2sf_subtree[asserting_cast<tsk_id_t>(sample_id)] = dag_node_id;
        }
        auto is_sample = [num_samples = _tree_sequence.num_samples()](tsk_id_t ts_node_id) {
            return ts_node_id < asserting_cast<tsk_id_t>(num_samples);
        };

        // TODO Rewrite this, once we have the tree_sequence iterator
        // TODO Add progress bar or progress report (for the large datasets)
        size_t const num_trees            = _tree_sequence.num_trees();
        size_t       tree_counter         = 0;
        size_t       report_every_n_trees = num_trees / 100;
        for (ts_tree.first(); ts_tree.is_valid(); ts_tree.next()) {
            tree_counter++;
            if (report_every_n_trees == 0 || tree_counter % report_every_n_trees == 0) {
                std::cerr << "Compressing tree " << tree_counter << "/" << num_trees << std::endl;
            }
            // Reset the mapper for this tree
            // TODO Specify for which nodes we want mapping instead of computing and storing it for all nodes
            ts_node2sf_subtree.reset();

            for (auto const ts_node_id: ts_tree.postorder()) {
                // Prefetching the children of the next node slows down the code in benchmarks.
                if (is_sample(ts_node_id)) {
                    // auto dag_subtree_id = succinct_subtree_id_factory.compute(ts_node_id);
                    // Samples are already mapped and added to this map above
                    auto dag_subtree_id = ts_node_to_dag_subtree_id_map[asserting_cast<size_t>(ts_node_id)];
                    KASSERT(
                        _dag_subtree_to_node_map.contains(dag_subtree_id),
                        "A leaf present in a later tree is missing from the first tree.",
                        tdt::assert::light
                    );
                    NodeId dag_node_id             = _dag_subtree_to_node_map.get(dag_subtree_id);
                    ts_node2sf_subtree[ts_node_id] = dag_node_id;
                } else { // Node is inner node
                    // Compute the subtree id of this inner node by hashing the subtree ids of it's children.
                    static_assert(std::is_same_v<decltype(XXH128_hash_t::low64), decltype(XXH128_hash_t::high64)>);
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
                        ts_node2sf_subtree[ts_node_id] = dag_node_id;
                    } else {
                        ts_node2sf_subtree[ts_node_id] = _dag_subtree_to_node_map[dag_subtree_id];
                    }
                }
            }

            // Process the mutations of this tree
            genomic_sequence_storage_factory.process_mutations(
                asserting_cast<TreeId>(ts_tree.tree_id()),
                ts_node2sf_subtree
            );
        }

        genomic_sequence_storage_factory.finalize();

        KASSERT(compressed_forest.roots().size() == _tree_sequence.num_trees());
        KASSERT(compressed_forest.num_roots() == _tree_sequence.num_trees());
        KASSERT(compressed_forest.num_leaves() == _tree_sequence.num_samples());

        // Set the number of nodes in the DAG so it does not have to be recomputed.
        compressed_forest.num_nodes(_dag_subtree_to_node_map.num_nodes());

        // As we build the tree edges by a postorder traversal on the tree, the from edges should be post-ordered, too.
        compressed_forest.postorder_edges().traversal_order(TraversalOrder::Postorder);

        return compressed_forest;
    }

private:
    SubtreeId2NodeMapper _dag_subtree_to_node_map;
    TSKitTreeSequence&   _tree_sequence;
};
