#include <cstddef>
#include <debug/vector>

#include <catch2/catch_test_macros.hpp>
#include <stdint.h>
#include <xxhash.h>

#include "kassert/kassert.hpp"
#include "sfkit/utils/xxhash.hpp"

TEST_CASE("XXH128_hash_t comparison", "[xxhash]") {
    XXH128_hash_t const hash1 = {0, 0};
    XXH128_hash_t const hash2 = {0, 0};
    XXH128_hash_t const hash3 = {1, 1};
    XXH128_hash_t const hash4 = {1, 2};
    XXH128_hash_t const hash5 = {1, 1};

    CHECK(hash1 == hash1);
    CHECK(hash2 == hash2);
    CHECK(hash3 == hash3);
    CHECK(hash4 == hash4);
    CHECK(hash5 == hash5);

    CHECK(hash2 == hash1);
    CHECK(hash1 == hash2);
    CHECK(hash3 == hash3);
    CHECK(hash5 == hash5);

    CHECK(hash1 != hash3);
    CHECK(hash3 != hash1);
    CHECK(hash1 != hash4);
    CHECK(hash4 != hash1);
    CHECK(hash4 != hash5);
    CHECK(hash5 != hash4);
}

TEST_CASE("xxhash128(POD)", "[xxhash]") {
    auto check_is_wrapper = [](auto const& value) {
        XXH64_hash_t const seed = 0;
        CHECK(xxhash128(value, seed) == XXH3_128bits_withSeed(&value, sizeof(value), seed));
    };

    check_is_wrapper(uint8_t(0));
    check_is_wrapper(int8_t(0));
    check_is_wrapper(uint16_t(0));
    check_is_wrapper(int16_t(0));
    check_is_wrapper(uint32_t(0));
    check_is_wrapper(int32_t(0));
    check_is_wrapper(uint64_t(0));
    check_is_wrapper(int64_t(0));
    check_is_wrapper(std::byte(0));
    check_is_wrapper(char(0));
    check_is_wrapper(static_cast<signed char>(0));
    check_is_wrapper(static_cast<unsigned char>(0));
}

TEST_CASE("xxhash128(Container)", "[xxhash]") {
    std::vector<std::byte> const     container1 = {std::byte(0), std::byte(1), std::byte(2), std::byte(3)};
    std::vector<std::byte> const     container2 = {std::byte(0), std::byte(1), std::byte(2), std::byte(3)};
    std::vector<XXH128_hash_t> const container3 = {{0, 0}, {1, 1}, {2, 2}, {3, 3}};
    std::vector<XXH128_hash_t> const container4 = {{0, 0}, {1, 1}, {2, 2}, {3, 3}};
    std::vector<XXH128_hash_t> const container5 = {{0, 1}, {1, 2}, {2, 3}, {4, 3}};
    std::vector<XXH128_hash_t> const container6 = {{0, 1}, {1, 2}, {2, 3}, {4, 3}};

    KASSERT(container1.size() == 4u);
    KASSERT(container2.size() == 4u);
    KASSERT(container3.size() == 4u);
    KASSERT(container4.size() == 4u);
    KASSERT(container5.size() == 4u);
    KASSERT(container6.size() == 4u);

    CHECK(xxhash128(container1) == xxhash128(container1));
    CHECK(xxhash128(container1) == xxhash128(container2));
    CHECK(xxhash128(container1) != xxhash128(container3));
    CHECK(xxhash128(container1) != xxhash128(container4));
    CHECK(xxhash128(container1) != xxhash128(container5));
    CHECK(xxhash128(container1) != xxhash128(container6));

    CHECK(xxhash128(container2) == xxhash128(container2));
    CHECK(xxhash128(container2) != xxhash128(container3));
    CHECK(xxhash128(container2) != xxhash128(container4));
    CHECK(xxhash128(container2) != xxhash128(container5));
    CHECK(xxhash128(container2) != xxhash128(container6));

    CHECK(xxhash128(container3) == xxhash128(container3));
    CHECK(xxhash128(container3) == xxhash128(container4));
    CHECK(xxhash128(container3) != xxhash128(container5));
    CHECK(xxhash128(container3) != xxhash128(container6));

    CHECK(xxhash128(container4) == xxhash128(container4));
    CHECK(xxhash128(container4) != xxhash128(container5));
    CHECK(xxhash128(container4) != xxhash128(container6));

    CHECK(xxhash128(container5) == xxhash128(container5));
    CHECK(xxhash128(container5) == xxhash128(container6));

    CHECK(xxhash128(container6) == xxhash128(container6));
}
