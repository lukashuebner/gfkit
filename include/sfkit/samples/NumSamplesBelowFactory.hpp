#pragma once

#include "sfkit/samples/NumSamplesBelow.hpp"
#include "sfkit/samples/NumSamplesBelowAccessor.hpp"
#include "sfkit/samples/SampleSet.hpp"

namespace sfkit::samples {

using sfkit::graph::EdgeListGraph;

class NumSamplesBelowFactory {
public:
    // TODO Use templates to avoid code duplication
    template <typename CompressedForest, typename BaseType = SampleId>
    static NumSamplesBelowAccessor<NumSamplesBelow<CompressedForest, 1, BaseType>>
    build(CompressedForest const& forest, SampleSet const& samples) {
        constexpr auto num_sample_sets = 1;
        using NumSamplesBelow          = NumSamplesBelow<CompressedForest, num_sample_sets, BaseType>;
        using SetOfSampleSets          = SetOfSampleSets<num_sample_sets>;

        auto const num_samples_below = std::make_shared<NumSamplesBelow>(forest, SetOfSampleSets{std::cref(samples)});

        return NumSamplesBelowAccessor(num_samples_below, 0);
    }

    template <typename CompressedForest, typename BaseType = SampleId>
    static auto build(CompressedForest const& forest, SampleSet const& samples_0, SampleSet const& samples_1) {
        constexpr auto num_sample_sets = 2;
        using NumSamplesBelow          = NumSamplesBelow<CompressedForest, num_sample_sets, BaseType>;
        using SetOfSampleSets          = SetOfSampleSets<num_sample_sets>;

        auto num_samples_below =
            std::make_shared<NumSamplesBelow>(forest, SetOfSampleSets{std::cref(samples_0), std::cref(samples_1)});

        return std::tuple(NumSamplesBelowAccessor(num_samples_below, 0), NumSamplesBelowAccessor(num_samples_below, 1));
    }

    template <typename CompressedForest, typename BaseType = SampleId>
    static auto build(
        CompressedForest const& forest,
        SampleSet const&        samples_0,
        SampleSet const&        samples_1,
        SampleSet const&        samples_2,
        SampleSet const&        samples_3
    ) {
        constexpr auto num_sample_sets = 4;
        using NumSamplesBelow          = NumSamplesBelow<CompressedForest, num_sample_sets, BaseType>;
        using SetOfSampleSets          = SetOfSampleSets<num_sample_sets>;

        auto num_samples_below = std::make_shared<NumSamplesBelow>(
            forest,
            SetOfSampleSets{std::cref(samples_0), std::cref(samples_1), std::cref(samples_2), std::cref(samples_3)}
        );

        return std::tuple(
            NumSamplesBelowAccessor(num_samples_below, 0),
            NumSamplesBelowAccessor(num_samples_below, 1),
            NumSamplesBelowAccessor(num_samples_below, 2),
            NumSamplesBelowAccessor(num_samples_below, 3)
        );
    }

private:
};

} // namespace sfkit::samples
