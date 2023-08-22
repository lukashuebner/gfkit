#pragma once

#include <string>
#include <typeinfo>

#include <kassert/kassert.hpp>
#include <xxhash.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/utils/xxhash.hpp"

#define TDT_XXHASH_SUBTREE_IDS

#ifdef TDT_XXHASH_SUBTREE_IDS
using SubtreeHash                           = XXH128_hash_t;
constexpr SubtreeHash SuccinctSubtreeIdZero = {0, 0};

inline SubtreeHash& operator^=(SubtreeHash& lhs, SubtreeHash const& rhs) {
    lhs.low64 ^= rhs.low64;
    lhs.high64 ^= rhs.high64;
    return lhs;
}

class SubtreeHasher {
public:
    SubtreeHasher(XXH64_hash_t const& seed = 42) : _seed(seed) {}

    template <typename T>
    requires requires(T const& t, XXH64_hash_t const& seed) {
        xxhash128(t, seed);
    }
    SubtreeHash compute(T const& data) {
        return xxhash128(data, _seed);
    }

    // Compute the subtree ID of this inner node by hashing the subtree IDs of its children.
    // template <typename Container>
    // XXH128_hash_t compute(Container const& children) {
    //     static_assert(std::is_same_v<decltype(XXH128_hash_t::low64), decltype(XXH128_hash_t::high64)>);
    //     SubtreeHash subtree_id = SuccinctSubtreeIdZero;

    //     // TODO Don't reserve memory for this. Save the current position in the edge list, add the target
    //     // nodes without the from edge and fix the from edge at the end of this block.
    //     std::vector<NodeId> children_dag_node_ids;

    //     for (auto&& child_id: children) {
    //         auto const child_dag_subtree_id = ts_node_to_subtree[asserting_cast<size_t>(child_ts_id)];

    //         children_dag_subtree_ids ^= child_dag_subtree_id;

    //         KASSERT(_subtree_to_sf_node.contains(child_dag_subtree_id));
    //         children_dag_node_ids.emplace_back(_subtree_to_sf_node[child_dag_subtree_id]);
    //     }
    // }

private:
    XXH64_hash_t const _seed;
};

// SubtreeHash already is a hash; this function thus is the identity.
template <>
struct std::hash<SubtreeHash> {
    size_t operator()(SubtreeHash const& subtree_id) const noexcept {
        return *(reinterpret_cast<size_t const*>(&subtree_id));
    }
};

#endif
