#pragma once

#include <array>

#include <kassert/kassert.hpp>
#include <xxhash.h>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/utils/concepts.hpp"

// TODO Test this for more data types (e.g. std::array)

template <typename Value>
requires std::has_unique_object_representations_v<Value>
    XXH128_hash_t xxhash128(Value const& value, XXH64_hash_t const& seed = 0) {
    static_assert(std::is_trivially_copyable_v<Value>);
    return XXH3_128bits_withSeed(&value, sizeof(Value), seed);
}

template <PlainStorageContainer Container>
requires std::has_unique_object_representations_v<typename Container::value_type>
    XXH128_hash_t xxhash128(Container const& container, XXH64_hash_t seed = 0) {
    size_t const num_bytes = sizeof(typename Container::value_type) * container.size();
    return XXH3_128bits_withSeed(container.data(), num_bytes, seed);
}

inline bool operator==(XXH128_hash_t const& lhs, XXH128_hash_t const& rhs) {
    return lhs.low64 == rhs.low64 && lhs.high64 == rhs.high64;
}

inline std::ostream& operator<<(std::ostream& os, XXH128_hash_t const& hash) {
    os << hash.low64 << hash.high64;
    return os;
}

template <typename Value>
requires std::has_unique_object_representations_v<Value>
    XXH64_hash_t xxhash64(Value const& value, XXH64_hash_t const& seed = 0) {
    static_assert(std::is_trivially_copyable_v<Value>);
    return XXH3_64bits_withSeed(&value, sizeof(Value), seed);
}

template <PlainStorageContainer Container>
requires std::has_unique_object_representations_v<typename Container::value_type>
    XXH64_hash_t xxhash64(Container const& container, XXH64_hash_t seed = 0) {
    size_t const num_bytes = sizeof(typename Container::value_type) * container.size();
    return XXH3_64bits_withSeed(container.data(), num_bytes, seed);
}
