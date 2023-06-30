#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include <tsl/hopscotch_map.h>
#pragma GCC diagnostic pop


#include "tdt/graph/common.hpp"
#include "tdt/load/subtree-id.hpp"

// TODO Rename files to use the same case as the class names
class SubtreeId2NodeMapper {
public:
    SubtreeId2NodeMapper() {
        // TODO Retest this with the new hash map
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
    // mutable google::dense_hash_map<SuccinctSubtreeId, NodeId> _subtree_to_node_map;
    mutable tsl::hopscotch_map<SuccinctSubtreeId, NodeId> _subtree_to_node_map;
    NodeId                                                _next_node_id = 0;
};