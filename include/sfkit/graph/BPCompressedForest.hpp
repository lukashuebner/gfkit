#pragma once

// #include <sparsehash/dense_hash_map>
#include <memory>
#include <unordered_set>

#include <kassert/kassert.hpp>
#include <sfkit/include-redirects/cereal.hpp>
#include <sfkit/include-redirects/sdsl.hpp>

// TODO include-what-you-use
#include "sfkit/assertion_levels.hpp"
#include "sfkit/checking_casts.hpp"
#include "sfkit/graph/AdjacencyArrayGraph.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/graph/balanced_parenthesis.hpp"
#include "sfkit/load/SubtreeHashToNodeMapper.hpp"
#include "sfkit/load/SubtreeHasher.hpp"
#include "sfkit/samples/NumSamplesBelow.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/tskit.hpp"
#include "sfkit/utils/concepts.hpp"

class BPCompressedForest {
public:
    static constexpr bool PARENS_OPEN  = bp::PARENS_OPEN;
    static constexpr bool PARENS_CLOSE = bp::PARENS_CLOSE;

    // TODO We don't move these vectors as we want to compress them. -> Add compression
    // TODO add const where possible
    BPCompressedForest(
        // TODO Use proper types in bp namespace
        sdsl::bit_vector const&                  is_reference,
        sdsl::bit_vector const&                  is_leaf,
        sdsl::bit_vector const&                  balanced_parenthesis,
        sdsl::int_vector<> const&                references,
        sdsl::int_vector<NodeId_bitwidth> const& leaves,
        NodeId const                             num_nodes,
        NodeId const                             num_leaves,
        TreeId const                             num_trees
    )
        : _is_reference(is_reference),
          _is_leaf(is_leaf),
          _balanced_parenthesis(balanced_parenthesis),
          _references(references),
          _leaves(leaves),
          _num_nodes(num_nodes),
          _num_leaves(num_leaves),
          _num_trees(num_trees) {
        sdsl::util::init_support(_is_reference_rank, &_is_reference);
        sdsl::util::init_support(_is_leaf_rank, &_is_leaf);
        sdsl::util::init_support(_balanced_parenthesis_rank, &_balanced_parenthesis);
    }

    [[nodiscard]] NodeId num_nodes() const {
        return _num_nodes;
    }

    [[nodiscard]] SampleSet all_samples() const {
        SampleSet sample_set(num_samples());
        for (NodeId sample_id = 0; sample_id < num_samples(); ++sample_id) {
            sample_set.add(sample_id);
        }
        return sample_set;
    }

    [[nodiscard]] SampleId num_samples() const {
        return _num_leaves;
    }

    [[nodiscard]] TreeId num_trees() const {
        return _num_trees;
    }

    [[nodiscard]] SampleId num_leaves() const {
        return _num_leaves;
    }

    [[nodiscard]] NodeId num_unique_subtrees() const {
        return _num_nodes;
    }

    // [[nodiscard]] EdgeId num_edges() const {
    //     return _dag_postorder_edges.num_edges();
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

    NodeId node_id(size_t const bp_idx) const {
        if (_is_leaf[bp_idx]) {
            auto const nth_leaf = _is_leaf_rank(bp_idx) >> 1;
            KASSERT(nth_leaf < num_leaves());
            std::cout << "bp_idx: " << bp_idx << " nth_leaf: " << nth_leaf
                      << " _leaves[nth_leaf]: " << _leaves[nth_leaf] << std::endl;
            return _leaves[nth_leaf];
        } else {
            auto const rank_in_bp              = _balanced_parenthesis_rank(bp_idx);
            auto const num_leaf_bits_in_prefix = _is_leaf_rank(bp_idx) >> 1;
            KASSERT(
                rank_in_bp >= num_leaf_bits_in_prefix >> 1,
                "Inconsistent state of rank support vectors.",
                sfkit::assert::light
            );
            std::cout << "bp_idx: " << bp_idx << " rank_in_bp: " << rank_in_bp
                      << " num_leaf_bits_in_prefix: " << num_leaf_bits_in_prefix << " num_leaves(): " << num_leaves()
                      << " node_id: " << rank_in_bp - num_leaf_bits_in_prefix + num_leaves() << std::endl;
            return asserting_cast<NodeId>(rank_in_bp - num_leaf_bits_in_prefix + num_leaves());
        }
    }

    NodeId node_id_ref(size_t const bp_idx) const {
        KASSERT((_is_reference_rank(bp_idx) & 1ul) == 0ul);
        // TODO document, why we're not dividing by 2
        auto const ref_rank         = _is_reference_rank(bp_idx);
        auto const idx_of_reference = _references[ref_rank];
        return node_id(idx_of_reference);
    }

    // This function is accessible mainly for unit-testing. It is not part of the public API.
    auto const& is_reference() const {
        return _is_reference;
    }

    auto const& balanced_parenthesis() const {
        return _balanced_parenthesis;
    }

    auto const& references() const {
        return _references;
    }

    auto const& leaves() const {
        return _leaves;
    }

    auto const& is_leaf() const {
        return _is_leaf;
    }

private:
    // TODO Use compressed bit-vectors
    sdsl::bit_vector        _is_reference;
    sdsl::rank_support_v5<> _is_reference_rank;
    sdsl::bit_vector        _is_leaf;
    sdsl::rank_support_v5<> _is_leaf_rank;
    sdsl::bit_vector        _balanced_parenthesis;
    // TODO Document why we're ranking the 0 pattern
    sdsl::rank_support_v5<0>          _balanced_parenthesis_rank;
    sdsl::int_vector<>                _references;
    sdsl::int_vector<NodeId_bitwidth> _leaves;
    NodeId                            _num_nodes;
    SampleId                          _num_leaves;
    TreeId                            _num_trees;
};
