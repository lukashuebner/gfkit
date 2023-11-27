#pragma once

#include "sfkit/bp/BPCompressedForest.hpp"
#include "sfkit/dag/DAGCompressedForest.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/samples/BPNumSamplesBelow.hpp"
#include "sfkit/samples/DAGNumSamplesBelow.hpp"
#include "sfkit/samples/NumSamplesBelowAccessor.hpp"
#include "sfkit/samples/SampleSet.hpp"

namespace sfkit::samples {

using sfkit::bp::BPCompressedForest;

class NumSamplesBelowFactory {
public:
    template <typename BaseType = SampleId>
    static NumSamplesBelowAccessor<DAGNumSamplesBelow<1, BaseType>>
    build(EdgeListGraph const& dag, SampleSet const& samples) {
        using SetOfSampleSets          = typename DAGNumSamplesBelow<1, BaseType>::SetOfSampleSets;
        constexpr auto num_sample_sets = 1;

        auto const num_samples_below =
            std::make_shared<DAGNumSamplesBelow<num_sample_sets, BaseType>>(dag, SetOfSampleSets{std::cref(samples)});

        return NumSamplesBelowAccessor(num_samples_below, 0);
    }

    template <typename BaseType = SampleId>
    static NumSamplesBelowAccessor<BPNumSamplesBelow<1, BaseType>>
    build(BPCompressedForest const& forest, SampleSet const& samples) {
        using SetOfSampleSets          = typename BPNumSamplesBelow<1, BaseType>::SetOfSampleSets;
        constexpr auto num_sample_sets = 1;

        auto const num_samples_below =
            std::make_shared<BPNumSamplesBelow<num_sample_sets, BaseType>>(forest, SetOfSampleSets{std::cref(samples)});

        return NumSamplesBelowAccessor(num_samples_below, 0);
    }

    template <typename BaseType = SampleId>
    static auto build(EdgeListGraph const& dag, SampleSet const& samples_0, SampleSet const& samples_1) {
        using SetOfSampleSets          = typename DAGNumSamplesBelow<2, BaseType>::SetOfSampleSets;
        constexpr auto num_sample_sets = 2;

        auto num_samples_below = std::make_shared<DAGNumSamplesBelow<num_sample_sets>>(
            dag,
            SetOfSampleSets{std::cref(samples_0), std::cref(samples_1)}
        );

        return std::tuple(NumSamplesBelowAccessor(num_samples_below, 0), NumSamplesBelowAccessor(num_samples_below, 1));
    }

    // TODO Use template magic to avoid code duplication
    template <typename BaseType = SampleId>
    static auto build(BPCompressedForest const& forest, SampleSet const& samples_0, SampleSet const& samples_1) {
        using SetOfSampleSets          = typename BPNumSamplesBelow<2, BaseType>::SetOfSampleSets;
        constexpr auto num_sample_sets = 2;

        auto num_samples_below = std::make_shared<BPNumSamplesBelow<num_sample_sets>>(
            forest,
            SetOfSampleSets{std::cref(samples_0), std::cref(samples_1)}
        );

        return std::tuple(NumSamplesBelowAccessor(num_samples_below, 0), NumSamplesBelowAccessor(num_samples_below, 1));
    }

    template <typename BaseType = SampleId>
    static auto build(
        EdgeListGraph const& dag,
        SampleSet const&     samples_0,
        SampleSet const&     samples_1,
        SampleSet const&     samples_2,
        SampleSet const&     samples_3
    ) {
        using SetOfSampleSets          = typename DAGNumSamplesBelow<4, BaseType>::SetOfSampleSets;
        constexpr auto num_sample_sets = 4;

        auto num_samples_below = std::make_shared<DAGNumSamplesBelow<num_sample_sets, BaseType>>(
            dag,
            SetOfSampleSets{std::cref(samples_0), std::cref(samples_1), std::cref(samples_2), std::cref(samples_3)}
        );

        return std::tuple(
            NumSamplesBelowAccessor(num_samples_below, 0),
            NumSamplesBelowAccessor(num_samples_below, 1),
            NumSamplesBelowAccessor(num_samples_below, 2),
            NumSamplesBelowAccessor(num_samples_below, 3)
        );
    }

    template <typename BaseType = SampleId>
    static auto build(
        BPCompressedForest const& forest,
        SampleSet const&          samples_0,
        SampleSet const&          samples_1,
        SampleSet const&          samples_2,
        SampleSet const&          samples_3
    ) {
        using SetOfSampleSets          = typename BPNumSamplesBelow<4, BaseType>::SetOfSampleSets;
        constexpr auto num_sample_sets = 4;

        auto num_samples_below = std::make_shared<BPNumSamplesBelow<num_sample_sets, BaseType>>(
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
};
} // namespace sfkit::samples
