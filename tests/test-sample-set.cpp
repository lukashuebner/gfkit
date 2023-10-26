
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_container_properties.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <tskit.h>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/samples/NumSamplesBelow.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch::Matchers;

TEST_CASE("SampleSet Basics", "[SampleSet]") {
    SampleSet samples{10};

    CHECK(samples.overall_num_samples() == 10);
    CHECK(samples.popcount() == 0);

    samples.add(0);
    samples.add(9);
    samples.add(5);

    CHECK(samples.overall_num_samples() == 10);
    CHECK(samples.popcount() == 3);

    // Check the operator[]
    CHECK(samples[0] == true);
    CHECK(samples[1] == false);
    CHECK(samples[2] == false);
    CHECK(samples[3] == false);
    CHECK(samples[4] == false);
    CHECK(samples[5] == true);
    CHECK(samples[6] == false);
    CHECK(samples[7] == false);
    CHECK(samples[8] == false);
    CHECK(samples[9] == true);

    // Check the iterator
    CHECK_THAT(samples, RangeEquals(std::vector<SampleId>{0, 5, 9}));
}

TEST_CASE("SampleSet::inverse()", "[SampleSet]") {
    SampleSet samples{10};

    CHECK(samples.overall_num_samples() == 10);
    CHECK(samples.inverse().overall_num_samples() == samples.overall_num_samples());
    CHECK(samples.inverse().popcount() == 10);

    samples.add(0);
    samples.add(9);
    samples.add(5);

    SampleSet inverse = samples.inverse();
    CHECK(inverse.overall_num_samples() == samples.overall_num_samples());
    CHECK(inverse.popcount() == 7);

    // Check the operator[]
    CHECK(inverse[0] == false);
    CHECK(inverse[1] == true);
    CHECK(inverse[2] == true);
    CHECK(inverse[3] == true);
    CHECK(inverse[4] == true);
    CHECK(inverse[5] == false);
    CHECK(inverse[6] == true);
    CHECK(inverse[7] == true);
    CHECK(inverse[8] == true);
    CHECK(inverse[9] == false);

    // Check the iterator
    CHECK_THAT(inverse, RangeEquals(std::vector<SampleId>{1, 2, 3, 4, 6, 7, 8}));
}

TEST_CASE("SampleSet::to_tsk_samples()", "[SampleSet]") {
    SampleSet samples{10};

    { // Empty SampleSet
        auto const tsk_samples = samples.to_tsk_samples();
        CHECK_THAT(tsk_samples, IsEmpty());
    }

    { // Example with some samples
        samples.add(0);
        samples.add(1);
        samples.add(3);
        CHECK_THAT(samples, RangeEquals(std::vector<SampleId>{0, 1, 3}));
    }

    { // Example with some more samples
        samples.add(0);
        samples.add(1);
        samples.add(3);
        samples.add(8);
        samples.add(9);
        CHECK_THAT(samples, RangeEquals(std::vector<SampleId>{0, 1, 3, 8, 9}));
    }

    { // Example with all samples
        for (uint32_t i = 0; i < 10; ++i) {
            samples.add(i);
        }
        CHECK_THAT(samples, RangeEquals(std::vector<SampleId>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}));
    }
}
