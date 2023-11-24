#include <debug/vector>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <stdint.h>
#include <stdlib.h>

#include "catch2/catch_assertion_info.hpp"
#include "sfkit/include-redirects/sdsl.hpp"
#include "sfkit/utils/BufferedSDSLBitVectorView.hpp"

using sfkit::utils::BufferedSDSLBitVectorView;

using namespace ::Catch;
using namespace ::Catch::Matchers;

TEST_CASE("BufferedSDSLBitVectorView Simple Example", "[Utils]") {
    sdsl::bit_vector vec(10, 0);

    vec[0] = 1;
    vec[1] = 1;
    vec[3] = 1;

    BufferedSDSLBitVectorView view(vec);

    CHECK_THAT(view, RangeEquals(std::vector<bool>{1, 1, 0, 1, 0, 0, 0, 0, 0, 0}));

    vec[0] = 0;
    vec[9] = 1;

    CHECK_THAT(view, RangeEquals(std::vector<bool>{0, 1, 0, 1, 0, 0, 0, 0, 0, 1}));
}

TEST_CASE("BufferedSDSLBitVectorView Various lengths", "[Utils]") {
    uint64_t const   num_bits = GENERATE(0ul, 1, 10, 100, 1000, 10000, 64, 128, 255, 256, 257);
    sdsl::bit_vector vec(num_bits, 0);

    std::vector<bool> expected;
    expected.reserve(num_bits);
    for (size_t idx = 0; idx < vec.size(); ++idx) {
        bool value = rand() % 2;
        vec[idx]   = value;
        expected.push_back(value);
    }

    BufferedSDSLBitVectorView view(vec);

    CHECK_THAT(view, RangeEquals(expected));
}
