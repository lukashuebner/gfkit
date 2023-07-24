#pragma once

#include <span>

#include <tskit/core.h>

#include "tdt/graph/common.hpp"
#include "tdt/sequence/sequence.hpp"

using MutationId = uint32_t;

class Mutation {
public:
    Mutation() = default;

    Mutation(
        SiteId site_id, TreeId tree_id, SubtreeId subtree_id, AllelicState state, AllelicState parent_state
    ) noexcept
        : _site_id(site_id),
          _parent_state(parent_state),
          _derived_state(state),
          _tree_id(tree_id),
          _subtree_id(subtree_id) {}

    [[nodiscard]] SiteId site_id() const {
        return _site_id;
    }

    [[nodiscard]] AllelicState allelic_state() const {
        return _derived_state;
    }

    [[nodiscard]] SubtreeId tree_id() const {
        return _tree_id;
    }

    [[nodiscard]] TreeId subtree_id() const {
        return _subtree_id;
    }

    template <class Archive>
    void serialize(Archive& archive) {
        archive(_site_id, _derived_state, _tree_id, _subtree_id, _parent_state);
    }

    bool operator==(Mutation const& other) const noexcept {
        return _site_id == other._site_id && _derived_state == other._derived_state && _tree_id == other._tree_id
               && _subtree_id == other._subtree_id && _parent_state == other._parent_state;
    }

    AllelicState parent_state() const {
        return _parent_state;
    }

private:
    // TODO Compress this representation. Through the use of indices we could get rid of the tree id
    // maybe even the site id, using a retrieval data structure? We also know the mapping site_id -> tree_id
    SiteId       _site_id;
    AllelicState _parent_state;
    AllelicState _derived_state;
    TreeId       _tree_id;
    SubtreeId    _subtree_id; // TODO Change to compressed forest node id
};

using MutationView = std::span<Mutation const>;
