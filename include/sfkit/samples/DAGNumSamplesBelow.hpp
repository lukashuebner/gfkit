#pragma once

#include <array>
#include <experimental/simd>
#include <vector>

#include <kassert/kassert.hpp>
#include <plf_stack.h>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/bp/BPCompressedForest.hpp"
#include "sfkit/dag/DAGCompressedForest.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/samples/NumSamplesBelow.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/samples/primitives.hpp"
#include "sfkit/utils/BufferedSDSLBitVectorView.hpp"

namespace sfkit::samples::internal {

using sfkit::dag::DAGCompressedForest;
using sfkit::graph::EdgeListGraph;
using sfkit::graph::NodeId;
namespace stdx = std::experimental;

template <size_t N, typename BaseType>
class DAGNumSamplesBelowImpl {
public:
    using SetOfSampleSets = std::array<std::reference_wrapper<SampleSet const>, N>;

    DAGNumSamplesBelowImpl(EdgeListGraph const& dag, SetOfSampleSets const& samples) : _dag(dag) {
        // Check inputs
        KASSERT(_dag.check_postorder(), "DAG edges are not post-ordered.", sfkit::assert::normal);
        KASSERT(_dag.num_nodes() >= _dag.num_leaves(), "DAG has less nodes than leaves.", sfkit::assert::light);
        for (auto sample_set: samples) {
            KASSERT(
                _dag.num_leaves() <= sample_set.get().overall_num_samples(),
                "Number of leaves in the DAG is greater than he number of overall samples representable in the subtree "
                "size object. (NOT the number of samples actually in the SampleSet)",
                sfkit::assert::light
            );
            KASSERT(
                sample_set.get().popcount() <= sample_set.get().overall_num_samples(),
                "Number of samples in sample set is greater than he number of overall samples representable in the "
                "subtree size object. (NOT the number of samples actually in the SampleSet)",
                sfkit::assert::light
            );
            // Will make redundant comparisons, but will only execute in DEBUG mode anyway
            for (auto other_sample_set: samples) {
                KASSERT(
                    sample_set.get().overall_num_samples() == other_sample_set.get().overall_num_samples(),
                    "Sample sets have different number of overall representable samples. (NOT the number of samples "
                    "actually "
                    "in the SampleSet)",
                    sfkit::assert::light
                );
            }
        }

        // Initialize _num_sample_in_sample_set
        auto num_samples_in_sample_set_it = _num_samples_in_sample_set.begin();
        auto samples_it                   = samples.begin();
        while (samples_it != samples.end()) {
            *num_samples_in_sample_set_it = samples_it->get().popcount();
            num_samples_in_sample_set_it++;
            samples_it++;
        }

        // Compute the subtree sizes
        _compute(samples);
    }

    [[nodiscard]] SampleId num_samples_below(NodeId node_id, SampleSetId sample_set_id) const {
        KASSERT(node_id < _subtree_sizes.size(), "Subtree ID out of bounds.", sfkit::assert::light);
        KASSERT(sample_set_id >= 0 && sample_set_id <= N, "Sample set ID invalid.", sfkit::assert::light);

        return _subtree_sizes[node_id][sample_set_id];
    }

    [[nodiscard]] SampleId operator()(NodeId node_id, SampleSetId sample_set_id) const {
        return this->num_samples_below(node_id, sample_set_id);
    }

    [[nodiscard]] SampleId num_nodes_in_dag() const {
        return asserting_cast<SampleId>(_dag.num_nodes());
    }

    [[nodiscard]] SampleId num_samples_in_dag() const {
        return asserting_cast<SampleId>(_dag.num_leaves());
    }

    [[nodiscard]] SampleId num_samples_in_sample_set(SampleSetId sample_set_id) const {
        KASSERT(sample_set_id >= 0 && sample_set_id <= N, "Sample set ID invalid.", sfkit::assert::light);
        return _num_samples_in_sample_set[sample_set_id];
    }

private:
    using simd_t = stdx::fixed_size_simd<BaseType, N>;
    EdgeListGraph const&    _dag; // As a post-order sorted edge list
    std::array<SampleId, N> _num_samples_in_sample_set;
    std::vector<simd_t>     _subtree_sizes;

    void _compute(SetOfSampleSets const& samples) {
        KASSERT(_subtree_sizes.size() == 0ul, "Subtree sizes already computed.", sfkit::assert::light);
        _subtree_sizes.resize(_dag.num_nodes(), 0);

        KASSERT(samples.size() == N);
        for (size_t sample_set_idx = 0; sample_set_idx < N; sample_set_idx++) {
            for (SampleId sample: samples[sample_set_idx].get()) {
                _subtree_sizes[sample][sample_set_idx] = 1;
            }
        }

        // Prefetch the first few edges
        constexpr auto num_edges_to_prefetch = 128;
        auto           prefetch_it           = _dag.begin();
        for (size_t i = 0; i < num_edges_to_prefetch && prefetch_it != _dag.end(); i++) {
            __builtin_prefetch(&_subtree_sizes[prefetch_it->from()], 0, 3);
            __builtin_prefetch(&_subtree_sizes[prefetch_it->to()], 0, 3);
            prefetch_it++;
        }

        // Compute the subtree sizes using a post-order traversal
        auto work_it = _dag.begin();
        while (prefetch_it != _dag.end()) {
            // Add all N bit counts bit-parallel
            _subtree_sizes[work_it->from()] += _subtree_sizes[work_it->to()];

            for (size_t sample_set_idx = 0; sample_set_idx < N; sample_set_idx++) {
                KASSERT(
                    (_subtree_sizes[work_it->from()][sample_set_idx]) <= _num_samples_in_sample_set[sample_set_idx],
                    "Number of samples below a node exceeds the number of samples in the tree sequence.",
                    sfkit::assert::light
                );
            }

            work_it++;
            __builtin_prefetch(&_subtree_sizes[prefetch_it->from()], 0, 3);
            __builtin_prefetch(&_subtree_sizes[prefetch_it->to()], 0, 3);
            prefetch_it++;
        }

        while (work_it != _dag.end()) {
            _subtree_sizes[work_it->from()] += _subtree_sizes[work_it->to()];

            for (size_t sample_set_idx = 0; sample_set_idx < N; sample_set_idx++) {
                KASSERT(
                    (_subtree_sizes[work_it->from()][sample_set_idx]) <= _num_samples_in_sample_set[sample_set_idx],
                    "Number of samples below a node exceeds the number of samples in the tree sequence.",
                    sfkit::assert::light
                );
            }

            work_it++;
        }
    }
};

} // namespace sfkit::samples::internal

namespace sfkit::samples {
using sfkit::dag::DAGCompressedForest;
using sfkit::graph::EdgeListGraph;

template <size_t N, typename BaseType>
class NumSamplesBelow<DAGCompressedForest, N, BaseType> {
public:
    using SetOfSampleSets = std::array<std::reference_wrapper<SampleSet const>, N>;

    NumSamplesBelow(DAGCompressedForest const& forest, SetOfSampleSets const& samples)
        : _impl(forest.postorder_edges(), samples) {}

    [[nodiscard]] SampleId num_samples_below(NodeId node_id, SampleSetId sample_set_id) const {
        return _impl.num_samples_below(node_id, sample_set_id);
    }

    [[nodiscard]] SampleId operator()(NodeId node_id, SampleSetId sample_set_id) const {
        return _impl(node_id, sample_set_id);
    }

    [[nodiscard]] SampleId num_nodes_in_dag() const {
        return _impl.num_nodes_in_dag();
    }

    [[nodiscard]] SampleId num_samples_in_dag() const {
        return _impl.num_samples_in_dag();
    }

    [[nodiscard]] SampleId num_samples_in_sample_set(SampleSetId sample_set_id) const {
        return _impl.num_samples_in_sample_set(sample_set_id);
    }

private:
    internal::DAGNumSamplesBelowImpl<N, BaseType> _impl;
};

template <size_t N, typename BaseType>
class NumSamplesBelow<EdgeListGraph, N, BaseType> {
public:
    using SetOfSampleSets = std::array<std::reference_wrapper<SampleSet const>, N>;

    NumSamplesBelow(EdgeListGraph const& dag, SetOfSampleSets const& samples) : _impl(dag, samples) {}

    [[nodiscard]] SampleId num_samples_below(NodeId node_id, SampleSetId sample_set_id) const {
        return _impl.num_samples_below(node_id, sample_set_id);
    }

    [[nodiscard]] SampleId operator()(NodeId node_id, SampleSetId sample_set_id) const {
        return _impl(node_id, sample_set_id);
    }

    [[nodiscard]] SampleId num_nodes_in_dag() const {
        return _impl.num_nodes_in_dag();
    }

    [[nodiscard]] SampleId num_samples_in_dag() const {
        return _impl.num_samples_in_dag();
    }

    [[nodiscard]] SampleId num_samples_in_sample_set(SampleSetId sample_set_id) const {
        return _impl.num_samples_in_sample_set(sample_set_id);
    }

private:
    internal::DAGNumSamplesBelowImpl<N, BaseType> _impl;
};

} // namespace sfkit::samples
