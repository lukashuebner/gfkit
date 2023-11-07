#pragma once

#include <sfkit/include-redirects/hopscotch_map.hpp>

#include "sfkit/graph/common.hpp"
#include "sfkit/load/SubtreeHasher.hpp"

class SubtreeHashToNodeMapper {
public:
    SubtreeHashToNodeMapper() {
        // TODO Retest this with the new hash map
        // Even if I know the exact size of the map, reserving the memory /degrades/ performance.
        // Hypothesis: Even more cache-misses in the beginning, when the map isn't fully filled yet.
        // _subtree_to_node_map.resize(...);
    }

    NodeId insert_node(SubtreeHash const& subtree_id) {
        KASSERT(!contains(subtree_id), "Subtree ID already exists in the map", sfkit::assert::light);
        _subtree_to_node_map[subtree_id] = _next_node_id;
        return _next_node_id++;
    }

    NodeId insert_or_update_node(SubtreeHash const& subtree_id) {
        _subtree_to_node_map[subtree_id] = _next_node_id;
        return _next_node_id++;
    }

    // Root nodes might have identical subtrees, thus the mapping is lo longer surjective. As roots are never
    // referred to, we do not need to store them but only assign them a node id.
    NodeId insert_root() {
        return _next_node_id++;
    }

    bool contains(SubtreeHash const& subtree_id) const {
        return _subtree_to_node_map.find(subtree_id) != _subtree_to_node_map.end();
    }

    auto find(SubtreeHash const& subtree_id) const {
        return _subtree_to_node_map.find(subtree_id);
    }

    auto end() const {
        return _subtree_to_node_map.end();
    }

    NodeId map(SubtreeHash const& subtree_id) const {
        KASSERT(contains(subtree_id), "Subtree ID does not exists in the map", sfkit::assert::light);
        return _subtree_to_node_map[subtree_id];
    }

    NodeId operator[](SubtreeHash const& subtree_id) const {
        return map(subtree_id);
    }

    NodeId num_nodes() const {
        return _next_node_id;
    }

private:
    using MapType = tsl::hopscotch_map<SubtreeHash, NodeId>;
    mutable MapType _subtree_to_node_map;
    NodeId          _next_node_id = 0;
};
