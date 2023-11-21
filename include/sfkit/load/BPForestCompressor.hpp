#pragma once

// #include <sparsehash/dense_hash_map>
#include <unordered_set>

#include <kassert/kassert.hpp>
#include <sfkit/include-redirects/hopscotch_map.hpp>
#include <sfkit/include-redirects/sdsl.hpp>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/checking_casts.hpp"
#include "sfkit/graph/AdjacencyArrayGraph.hpp"
#include "sfkit/graph/BPCompressedForest.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/graph/balanced_parenthesis.hpp"
#include "sfkit/load/SubtreeHashToNodeMapper.hpp"
#include "sfkit/load/SubtreeHasher.hpp"
#include "sfkit/load/TsToSfNodeMapper.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"
#include "sfkit/sequence/GenomicSequenceFactory.hpp"
#include "sfkit/tskit.hpp"
#include "sfkit/utils/concepts.hpp"

// TODO Separate this class into the construction and the storage
class BPForestCompressor {
public:
    BPForestCompressor(TSKitTreeSequence& tree_sequence)
        : _num_trees(tree_sequence.num_trees()),
          _ts_tree(tree_sequence) {
        if (!tree_sequence.sample_ids_are_consecutive()) {
            throw std::runtime_error("Sample IDs of the tree sequence are not consecutive.");
        }
        _ts_node_to_subtree.resize(_ts_tree.max_node_id());
        _register_samples(tree_sequence);
        KASSERT(_subtree_to_sf_node.num_nodes() == tree_sequence.num_samples());
    }

    template <typename GenomicSequenceFactoryT>
    BPCompressedForest compress(GenomicSequenceFactoryT& genomic_sequence_storage_factory) {
        _num_samples = 0;

        // TODO Rewrite this, once we have the tree_sequence iterator
        for (_ts_tree.first(); _ts_tree.is_valid(); _ts_tree.next()) {
            auto const eulertour = _ts_tree.eulertour();
            auto       node_it   = eulertour.begin();
            while (node_it != eulertour.end()) {
                auto const ts_node_id = node_it.node_id();
                if (node_it.is_sample()) {
                    if (_ts_tree.tree_id() == 0) {
                        _add_sample(asserting_cast<SampleId>(ts_node_id));
                    } else {
                        auto const subtree_id   = _ts_node_to_subtree[asserting_cast<size_t>(ts_node_id)];
                        auto const reference_it = _subtrees.find(subtree_id);
                        _refer_to(reference_it);
                    }
                } else if (node_it.first_visit()) {
                    // Add the subtree to the BP forest.
                    _open_subtree();
                } else {
                    KASSERT(node_it.second_visit());
                    // Compute the subtree ID of this inner node by hashing the XOR of the subtree IDs of its children.
                    _subtree_hash_factory.reset();

                    for (auto child_ts_id: _ts_tree.children(ts_node_id)) {
                        SubtreeHash const childs_subtree_id = _ts_node_to_subtree[asserting_cast<size_t>(child_ts_id)];

                        _subtree_hash_factory.append_child(childs_subtree_id);
                    }
                    auto const subtree_id = _subtree_hash_factory.hash();

                    // Did we already encounter this subtree and can refer to its encoding?
                    auto const reference_it = _subtrees.find(subtree_id);
                    if (reference_it != _subtrees.end()) {
                        // The referenced node (== subtree) has already been added to the BP, reference it.
                        _rollback_subtree();
                        _refer_to(reference_it);
                    } else {
                        // This is no subtree not encoded before, close and store if for future reference.
                        _close_and_commit_subtree(subtree_id);
                    }

                    // Multiple inner nodes in the tskit tree sequence might describe the same subtree. We recognize
                    // this and refer back to the already encoded node. We thus have to store the mapping ts node ->
                    // subtree (and thus ts node) if we encode a new subtree AND if we refer to an existing subtree.
                    _ts_node_to_subtree[asserting_cast<size_t>(ts_node_id)] = subtree_id;
                }
                ++node_it;
            }

            // Process the mutations of this tree
            genomic_sequence_storage_factory.process_mutations(
                asserting_cast<TreeId>(_ts_tree.tree_id()),
                TsToSfNodeMapper(_ts_node_to_subtree, _subtree_to_sf_node)
            );
        }

        genomic_sequence_storage_factory.finalize();

        _is_reference.shrink_to_fit();
        _is_leaf.shrink_to_fit();
        _balanced_parenthesis.shrink_to_fit();
        _references.shrink_to_fit();
        _leaves.shrink_to_fit();

        KASSERT(_subtrees.size() == _subtree_to_sf_node.num_nodes());
        KASSERT(_is_leaf.index() == _balanced_parenthesis.index());
        KASSERT(_is_reference.index() == _balanced_parenthesis.index());
        // TODO KASSERT popcounts of is_leaf == num entries in _leaves and _is_reference == num entries in _references
        return BPCompressedForest(
            _is_reference.underlying(),
            _is_leaf.underlying(),
            _balanced_parenthesis.underlying(),
            _references.underlying(),
            _leaves.underlying(),
            asserting_cast<NodeId>(_subtrees.size()),
            _num_samples, // TODO Do I need this if I have the leaves vector?
            _num_trees
        );
    }

private:
    struct Reference {
        size_t start;
        size_t length;
        NodeId node_id;
    };

    struct SubtreeStart {
        size_t bp_start;
        size_t ref_start;
        size_t leaves_start;
    };

    // TODO Move to another class
    template <typename SDSLContainer>
    class DynamicSDSLContainer {
    public:
        using size_type  = typename SDSLContainer::size_type;
        using value_type = typename SDSLContainer::value_type;

        DynamicSDSLContainer() : _container(initial_size) {}

        void push_back(value_type const& value) {
            if (_idx >= _container.size()) [[unlikely]] {
                _grow();
            }
            _container[_idx] = value;
            ++_idx;
        }

        void rollback_to(size_type const idx) {
            _idx = idx;
        }

        void shrink_to_fit() {
            _container.resize(_idx);
        }

        SDSLContainer const& underlying() const {
            return _container;
        }

        [[nodiscard]] size_type index() const {
            return _idx;
        }

    private:
        static constexpr double    grow_factor  = 1.95;
        static constexpr size_type initial_size = 1024;

        SDSLContainer _container;
        size_type     _idx = 0;

        void _grow() {
            _container.resize(static_cast<size_type>(grow_factor * static_cast<double>(_container.size())));
        }
    };

    using IsReference           = DynamicSDSLContainer<sdsl::bit_vector>;
    using IsLeaf                = DynamicSDSLContainer<sdsl::bit_vector>;
    using IsBalancedParenthesis = DynamicSDSLContainer<sdsl::bit_vector>;
    using References            = DynamicSDSLContainer<sdsl::int_vector<NodeId_bitwidth> >;
    using Leaves                = DynamicSDSLContainer<sdsl::int_vector<NodeId_bitwidth> >;
    using SubtreeStarts         = std::vector<SubtreeStart>;
    using Subtrees              = tsl::hopscotch_map<SubtreeHash, Reference>;

    SampleId                 _num_samples = 0;
    TreeId                   _num_trees   = 0;
    TSKitTree                _ts_tree;
    std::vector<SubtreeHash> _ts_node_to_subtree;
    SubtreeHashToNodeMapper  _subtree_to_sf_node;
    SubtreeHasher            _subtree_hash_factory;
    IsReference              _is_reference;
    IsLeaf                   _is_leaf;
    IsBalancedParenthesis    _balanced_parenthesis;
    References               _references;
    Leaves                   _leaves;
    SubtreeStarts            _subtree_starts;
    Subtrees                 _subtrees;

    // The sample ids are consecutive: 0 ... num_samples - 1
    inline bool is_sample(tsk_id_t ts_node_id) const {
        return ts_node_id < asserting_cast<tsk_id_t>(_num_samples);
    }

    void _add_sample(SampleId const sample_id) {
        _open_subtree(sample_id);
        SubtreeHash const subtree_id = _subtree_hash_factory.hash_sample(sample_id);
        _close_and_commit_subtree(subtree_id, true);
        ++_num_samples;
    }

    void _open_subtree(SampleId sample_id = INVALID_NODE_ID) {
        KASSERT(_balanced_parenthesis.index() == _is_reference.index());
        KASSERT(_balanced_parenthesis.index() == _is_leaf.index());
        _subtree_starts.emplace_back(_balanced_parenthesis.index(), _references.index(), _leaves.index());
        _is_reference.push_back(false);
        _balanced_parenthesis.push_back(bp::PARENS_OPEN);
        if (sample_id != INVALID_NODE_ID) {
            _is_leaf.push_back(true);
            _leaves.push_back(sample_id);
        } else {
            _is_leaf.push_back(false);
        }
        KASSERT(_balanced_parenthesis.index() == _is_reference.index());
        KASSERT(_balanced_parenthesis.index() == _is_leaf.index());
    }

    void _close_and_commit_subtree(SubtreeHash subtree_id, bool is_leaf = false) {
        KASSERT(_balanced_parenthesis.index() == _is_reference.index());
        KASSERT(_balanced_parenthesis.index() == _is_leaf.index());
        _is_reference.push_back(false);
        _balanced_parenthesis.push_back(bp::PARENS_CLOSE);
        _is_leaf.push_back(is_leaf);

        // TODO This should not be independent of the node_id calculation in BPCompressedForest
        NodeId node_id = 0;
        if (is_leaf) { // Leaves are already registered
            node_id = _subtree_to_sf_node[subtree_id];
        } else { // Inner nodes need to be registered and get assigned an ID
            node_id = _subtree_to_sf_node.insert_node(subtree_id);
        }
        KASSERT(_balanced_parenthesis.index() == _is_reference.index());
        KASSERT(_balanced_parenthesis.index() == _is_leaf.index());

        // Add the reference to the subtree that begins at the open parenthesis.
        auto const start  = _subtree_starts.back().bp_start;
        auto const length = _balanced_parenthesis.index() - start;
        _subtrees.emplace(subtree_id, Reference{start, length, node_id});
    }

    void _rollback_subtree() {
        KASSERT(_balanced_parenthesis.index() == _is_reference.index());
        KASSERT(_balanced_parenthesis.index() == _is_leaf.index());
        auto const& subtree_start = _subtree_starts.back();
        _is_reference.rollback_to(subtree_start.bp_start);
        _is_leaf.rollback_to(subtree_start.bp_start);
        _balanced_parenthesis.rollback_to(subtree_start.bp_start);
        _leaves.rollback_to(subtree_start.leaves_start);
        _references.rollback_to(subtree_start.ref_start);
        _subtree_starts.pop_back();
        KASSERT(_balanced_parenthesis.index() == _is_reference.index());
        KASSERT(_balanced_parenthesis.index() == _is_leaf.index());
    }

    void _refer_to(decltype(_subtrees)::const_iterator reference_it) {
        KASSERT(_balanced_parenthesis.index() == _is_reference.index());
        KASSERT(_balanced_parenthesis.index() == _is_leaf.index());
        _is_reference.push_back(true);
        _is_reference.push_back(true);
        _balanced_parenthesis.push_back(bp::PARENS_OPEN);
        _balanced_parenthesis.push_back(bp::PARENS_CLOSE);

        _is_leaf.push_back(false);
        _is_leaf.push_back(false);

        KASSERT(
            reference_it != _subtrees.end(),
            "Trying to refer to subtree which has not been stored yet.",
            sfkit::assert::light
        );
        // TODO Do we need these values at any point?
        // _references.push_back(reference_it->second.start);
        // _references.push_back(reference_it->second.length);
        _references.push_back(reference_it->second.node_id);
        // TODO Re-use the tree topology with different leaf names.
        KASSERT(_balanced_parenthesis.index() == _is_reference.index());
        KASSERT(_balanced_parenthesis.index() == _is_leaf.index());
    }

    // Assign subtree-ids to the samples. Do this first, before adding any other nodes, ensuring that samples always map
    // to the same subtree-ids.
    void _register_samples(TSKitTreeSequence& tree_sequence) {
        KASSERT(tree_sequence.sample_ids_are_consecutive(), "Sample IDs are not consecutive.");
        for (SampleId sample_id = 0; sample_id < tree_sequence.num_samples(); sample_id++) {
            KASSERT(tree_sequence.is_sample(asserting_cast<tsk_id_t>(sample_id)));

            // Compute the subtree ID of this sample (leaf) node by hashing its label.
            auto subtree_hash = _subtree_hash_factory.hash_sample(sample_id);

            // Cache the subtree ID for this ts node
            KASSERT(asserting_cast<size_t>(sample_id) < _ts_node_to_subtree.size());
            _ts_node_to_subtree[asserting_cast<size_t>(sample_id)] = subtree_hash;

            // Map the subtree ID to the corresponding node ID in the DAG.
            _subtree_to_sf_node.insert_node(subtree_hash);
        }
    }
};
