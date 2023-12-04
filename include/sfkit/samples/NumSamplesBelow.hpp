#pragma once

#include <cstddef>
#include <cstdint>

#include "sfkit/graph/primitives.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/samples/primitives.hpp"

namespace sfkit::samples {

using sfkit::graph::NodeId;
using sfkit::samples::SampleSetId;

template <typename T>
concept NumSamplesBelowC = requires(T t) {
    { t(NodeId{}, SampleSetId{}) } -> std::convertible_to<SampleId>;
    { t.num_nodes_in_dag() } -> std::convertible_to<SampleId>;
    { t.num_samples_in_dag() } -> std::convertible_to<SampleId>;
    { t.num_samples_in_sample_set() } -> std::convertible_to<SampleId>;
};

template <typename CompressedForest, size_t N = 1, typename BaseType = SampleId>
class NumSamplesBelow {
public:
    SampleId operator()(NodeId node_id, SampleSetId sample_set_id);
    SampleId num_nodes_in_dag();
    SampleId num_samples_in_dag();
    SampleId num_samples_in_sample_set();
};

} // namespace sfkit::samples
