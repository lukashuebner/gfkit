#pragma once

#include <cstddef>
#include <memory>
#include <unordered_set>

#include <fmt/core.h>
#include <fmt/format.h>
#include <kassert/kassert.hpp>
#include <sfkit/include-redirects/cereal.hpp>
#include <sfkit/include-redirects/sdsl.hpp>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/bp/Parens.hpp"
#include "sfkit/graph/AdjacencyArrayGraph.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/graph/SubtreeHashToNodeMapper.hpp"
#include "sfkit/graph/SubtreeHasher.hpp"
#include "sfkit/graph/primitives.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/tskit/tskit.hpp"
#include "sfkit/utils/checking_casts.hpp"
#include "sfkit/utils/concepts.hpp"

namespace sfkit::bp {
using namespace sfkit::graph;
using sfkit::samples::SampleId;  // TODO Remove this dependency
using sfkit::samples::SampleSet; // TODO Remove this dependency

class BPCompressedForest {
public:
    static constexpr bool PARENS_OPEN  = bp::PARENS_OPEN;
    static constexpr bool PARENS_CLOSE = bp::PARENS_CLOSE;

    // TODO Cereal has a proper way to construct a class which has no default constructor
    BPCompressedForest() = default;

    // TODO We don't move these vectors as we want to compress them. -> Add compression
    // TODO add const where possible
    BPCompressedForest(
        // TODO Use proper types in bp namespace
        sdsl::bit_vector const&                  is_reference,
        sdsl::bit_vector const&                  is_leaf,
        sdsl::bit_vector const&                  balanced_parenthesis,
        sdsl::int_vector<NodeId_bitwidth> const& references,
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

    [[nodiscard]] NodeId num_backrefs() const {
        return asserting_cast<NodeId>(_is_reference_rank.rank(_is_reference_rank.size() - 1));
    }

    bool is_leaf(size_t const bp_idx) const {
        return _is_leaf[bp_idx];
    }

    SampleId leaf_idx_to_id(SampleId const leaf_id) const {
        KASSERT(leaf_id < num_leaves());
        return _leaves[asserting_cast<size_t>(leaf_id)];
    }

    // TODO Add documentation, especially about runtime
    NodeId node_id(size_t const bp_idx) const {
        if (_is_leaf[bp_idx]) {
            SampleId const nth_leaf = asserting_cast<SampleId>(_is_leaf_rank(bp_idx) >> 1);
            KASSERT(nth_leaf < num_leaves());
            return leaf_idx_to_id(nth_leaf);
        } else {
            auto const rank_in_bp              = _balanced_parenthesis_rank(bp_idx);
            auto const num_leaf_bits_in_prefix = _is_leaf_rank(bp_idx) >> 1;
            auto const num_ref_bits_in_prefix  = _is_reference_rank(bp_idx) >> 1;
            KASSERT(
                rank_in_bp >= num_leaf_bits_in_prefix + num_ref_bits_in_prefix,
                "Inconsistent state of rank support vectors.",
                sfkit::assert::light
            );
            auto const node_id =
                asserting_cast<NodeId>(rank_in_bp - num_leaf_bits_in_prefix - num_ref_bits_in_prefix + num_leaves());
            KASSERT(node_id < num_nodes());
            return node_id;
        }
    }

    NodeId node_id_ref(size_t const bp_idx) const {
        KASSERT((_is_reference_rank(bp_idx) & 1ul) == 0ul);
        // TODO document, why we're dividing by 2 and multiplying by 3 -> Abstract away to ReferenceStorage
        // TODO Count references to get rid of the rank
        // TODO Do we need the other information about the references? -> different vectors (easier indexing)?
        auto const ref_rank = _is_reference_rank(bp_idx) >> 1;
        return node_id_ref_by_rank(ref_rank);
    }

    NodeId node_id_ref_by_rank(size_t const ref_rank) const {
        return _references[ref_rank];
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

    size_t bp_size_bit() const {
        return balanced_parenthesis().bit_size();
    }

    size_t is_ref_size_bit() const {
        return is_reference().bit_size();
    }

    size_t is_leaf_size_bit() const {
        return is_leaf().bit_size();
    }

    size_t references_size_bit() const {
        return references().bit_size();
    }

    size_t leaves_size_bit() const {
        return leaves().bit_size();
    }

    // TODO move this to free functions?
    void save(std::ostream& os) const {
        _is_reference.serialize(os);
        //_is_reference_rank.serialize(os);
        _is_leaf.serialize(os);
        //_is_leaf_rank.serialize(os);
        _balanced_parenthesis.serialize(os);
        //_balanced_parenthesis_rank.serialize(os);
        _references.serialize(os);
        _leaves.serialize(os);

        os.write(reinterpret_cast<char const*>(&_num_nodes), sizeof(_num_nodes));
        os.write(reinterpret_cast<char const*>(&_num_leaves), sizeof(_num_leaves));
        os.write(reinterpret_cast<char const*>(&_num_trees), sizeof(_num_trees));
    }

    void load(std::istream& is) {
        _is_reference.load((is));
        // _is_reference_rank.load((is));
        _is_leaf.load((is));
        // _is_leaf_rank.load((is));
        _balanced_parenthesis.load(is);
        // _balanced_parenthesis_rank.load(is);
        _references.load((is));
        _leaves.load((is));

        is.read(reinterpret_cast<char*>(&_num_nodes), sizeof(_num_nodes));
        is.read(reinterpret_cast<char*>(&_num_leaves), sizeof(_num_leaves));
        is.read(reinterpret_cast<char*>(&_num_trees), sizeof(_num_trees));

        // TODO Can't these be serialized and deserialized?
        sdsl::util::init_support(_is_reference_rank, &_is_reference);
        sdsl::util::init_support(_is_leaf_rank, &_is_leaf);
        sdsl::util::init_support(_balanced_parenthesis_rank, &_balanced_parenthesis);
    }

    [[nodiscard]] bool operator==(BPCompressedForest const& other) {
        return _is_reference == other._is_reference && _is_leaf == other._is_leaf
               && _balanced_parenthesis == other._balanced_parenthesis && _references == other._references
               && _leaves == other._leaves && _num_nodes == other._num_nodes && _num_leaves == other._num_leaves
               && _num_trees == other._num_trees;
    }

private:
    // TODO Which of these vectors are actually needed?
    // TODO Use compressed bit-vectors
    sdsl::bit_vector        _is_reference;
    sdsl::rank_support_v5<> _is_reference_rank;
    sdsl::bit_vector        _is_leaf;
    sdsl::rank_support_v5<> _is_leaf_rank;
    sdsl::bit_vector        _balanced_parenthesis;
    // TODO Document why we're ranking the 0 pattern
    sdsl::rank_support_v5<0>          _balanced_parenthesis_rank;
    sdsl::int_vector<NodeId_bitwidth> _references;
    sdsl::int_vector<NodeId_bitwidth> _leaves;
    NodeId                            _num_nodes;
    SampleId                          _num_leaves;
    TreeId                            _num_trees;
};

} // namespace sfkit::bp
