#pragma once

#include <string>
#include <typeinfo>

#include <kassert/kassert.hpp>
#include <xxhash.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/utils/sha256.hpp"
#include "tdt/utils/xxhash.hpp"

#define TDT_XXHASH_SUBTREE_IDS

#ifdef TDT_XXHASH_SUBTREE_IDS
using SuccinctSubtreeId = XXH128_hash_t;
constexpr SuccinctSubtreeId SuccinctSubtreeIdZero = {0, 0};

inline SuccinctSubtreeId& operator^= (SuccinctSubtreeId& lhs, SuccinctSubtreeId const& rhs) {
    lhs.low64 ^= rhs.low64;
    lhs.high64 ^= rhs.high64;
    return lhs;
}

class SuccinctSubtreeIdFactory {
public:
    SuccinctSubtreeIdFactory(XXH64_hash_t const& seed = 42) : _seed(seed) {}

    template <typename T>
    requires requires(T const& t, XXH64_hash_t const& seed) {
        xxhash128(t, seed);
    }
    XXH128_hash_t compute(T const& data) {
        return xxhash128(data, _seed);
    }

private:
    XXH64_hash_t const _seed;
};

template <>
struct std::hash<SuccinctSubtreeId> {
    size_t operator()(SuccinctSubtreeId const& subtree_id) const noexcept {
        return *(reinterpret_cast<size_t const*>(&subtree_id));
    }
};

#endif
