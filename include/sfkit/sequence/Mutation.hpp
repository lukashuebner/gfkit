#pragma once

#include <span>

#include <tskit/core.h>

#include "sfkit/graph/common.hpp"
#include "sfkit/sequence/Sequence.hpp"

using MutationId = uint32_t;

class Mutation {
public:
    Mutation() = default;

    Mutation(SiteId site_id, TreeId tree_id, NodeId node_id, AllelicState state, AllelicState parent_state) noexcept
        : _site_id(site_id),
          _parent_state(parent_state),
          _derived_state(state),
          _tree_id(tree_id),
          _node_id(node_id) {}

    [[nodiscard]] SiteId site_id() const {
        return _site_id;
    }

    [[nodiscard]] AllelicState allelic_state() const {
        return _derived_state;
    }

    [[nodiscard]] NodeId node_id() const {
        return _node_id;
    }

    template <class Archive>
    void serialize(Archive& archive) {
        archive(_site_id, _derived_state, _tree_id, _node_id, _parent_state);
    }

    bool operator==(Mutation const& other) const noexcept {
        return _site_id == other._site_id && _derived_state == other._derived_state && _tree_id == other._tree_id
               && _node_id == other._node_id && _parent_state == other._parent_state;
    }

    AllelicState parent_state() const {
        return _parent_state;
    }

private:
    // TODO Compress this representation. Through the use of indices we could get rid of the tree id
    // maybe even the site id, using a retrieval data structure? We also know the mapping site_id -> tree_id
    // TODO Maybe we don't need the parent and derived state at least for biallelic sites (either filter the mutation or
    // it toggles the state)
    SiteId       _site_id;
    AllelicState _parent_state;
    AllelicState _derived_state;
    TreeId       _tree_id;
    // TODO Directly encode the formula to compute the subtree size here. E.g. a bitmask of which samples ares below
    // this node. (xor with sample set and then do a popcount)
    NodeId _node_id;
};

using MutationView = std::span<Mutation const>;
