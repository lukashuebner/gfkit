#pragma once

#include <concepts>
#include <memory>

#include "sfkit/graph/primitives.hpp"
#include "sfkit/samples/BPNumSamplesBelow.hpp"
#include "sfkit/samples/DAGNumSamplesBelow.hpp"
#include "sfkit/samples/NumSamplesBelow.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/samples/primitives.hpp"

namespace sfkit::samples {

using sfkit::graph::NodeId;

template <typename T>
concept NumSamplesBelowAccessorC = requires(T t) {
    { t.num_samples_below(NodeId{}) } -> std::convertible_to<SampleId>;
    { t(NodeId{}) } -> std::convertible_to<SampleId>;
    { t[NodeId{}] } -> std::convertible_to<SampleId>;
    { t.num_nodes_in_dag() } -> std::convertible_to<SampleId>;
    { t.num_samples_in_dag() } -> std::convertible_to<SampleId>;
    { t.num_samples_in_sample_set() } -> std::convertible_to<SampleId>;
};

template <typename NumSamplesBelow>
class NumSamplesBelowAccessor {
public:
    NumSamplesBelowAccessor(std::shared_ptr<NumSamplesBelow> const& num_samples_below, SampleSetId sample_set_id)
        : _num_samples_below(num_samples_below),
          _sample_set_id(sample_set_id) {}

    [[nodiscard]] SampleId num_samples_below(NodeId node_id) const {
        return _num_samples_below->operator()(node_id, _sample_set_id);
    }

    [[nodiscard]] SampleId operator()(NodeId node_id) const {
        return this->num_samples_below(node_id);
    }

    [[nodiscard]] SampleId operator[](NodeId node_id) const {
        return this->num_samples_below(node_id);
    }

    [[nodiscard]] SampleId num_nodes_in_dag() const {
        return _num_samples_below->num_nodes_in_dag();
    }

    [[nodiscard]] SampleId num_samples_in_dag() const {
        return _num_samples_below->num_samples_in_dag();
    }

    [[nodiscard]] SampleId num_samples_in_sample_set() const {
        return _num_samples_below->num_samples_in_sample_set(_sample_set_id);
    }

private:
    std::shared_ptr<NumSamplesBelow> _num_samples_below;
    SampleSetId                      _sample_set_id;
};

} // namespace sfkit::samples
