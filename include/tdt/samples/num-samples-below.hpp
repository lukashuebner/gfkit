#pragma once

#include <array>
#include <experimental/simd>
#include <vector>

#include <kassert/kassert.hpp>

#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/samples/sample-set.hpp"
// TODO Unify naming of files
#include "tdt/assertion_levels.hpp"
#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/samples/sample-set.hpp"

namespace stdx = std::experimental;

template <typename T>
concept NumSamplesBelowAccessorC = requires(T t) {
    { t.num_samples_below(NodeId{}) } -> std::convertible_to<SampleId>;
    { t(NodeId{}) } -> std::convertible_to<SampleId>;
    { t[NodeId{}] } -> std::convertible_to<SampleId>;
    { t.num_nodes_in_dag() } -> std::convertible_to<SampleId>;
    { t.num_samples_in_dag() } -> std::convertible_to<SampleId>;
    { t.num_samples_in_sample_set() } -> std::convertible_to<SampleId>;
};

using SampleSetId = uint8_t;

template <size_t N = 1, typename BaseType = SampleId>
class NumSamplesBelow {
public:
    using SetOfSampleSets = std::array<std::reference_wrapper<SampleSet const>, N>;

    NumSamplesBelow(EdgeListGraph const& dag_postorder_edges, SetOfSampleSets const& samples)
        : _dag(dag_postorder_edges) {
        // Check inputs
        KASSERT(_dag.check_postorder(), "DAG edges are not post-ordered.", tdt::assert::normal);
        KASSERT(_dag.num_nodes() >= _dag.num_leaves(), "DAG has less nodes than leaves.", tdt::assert::light);
        for (auto sample_set: samples) {
            KASSERT(
                _dag.num_leaves() <= sample_set.get().overall_num_samples(),
                "Number of leaves in the DAG is greater than he number of overall samples representable in the subtree "
                "size object. (NOT the number of samples actually in the SampleSet)",
                tdt::assert::light
            );
            KASSERT(
                sample_set.get().popcount() <= sample_set.get().overall_num_samples(),
                "Number of samples in sample set is greater than he number of overall samples representable in the "
                "subtree size object. (NOT the number of samples actually in the SampleSet)",
                tdt::assert::light
            );
            // Will make redundant comprisons, but will only execute in DEBUG mode anyway
            for (auto other_sample_set: samples) {
                KASSERT(
                    sample_set.get().overall_num_samples() == other_sample_set.get().overall_num_samples(),
                    "Sample sets have different number of overall representable samples. (NOT the number of samples "
                    "actually "
                    "in the SampleSet)",
                    tdt::assert::light
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
        KASSERT(node_id < _subtree_sizes.size(), "Subtree ID out of bounds.", tdt::assert::light);
        KASSERT(sample_set_id >= 0 && sample_set_id <= N, "Sample set ID invalid.", tdt::assert::light);

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
        KASSERT(sample_set_id >= 0 && sample_set_id <= N, "Sample set ID invalid.", tdt::assert::light);
        return _num_samples_in_sample_set[sample_set_id];
    }

private:
    using simd_t = stdx::fixed_size_simd<BaseType, N>;
    EdgeListGraph const&    _dag; // As a post-order sorted edge list
    std::array<SampleId, N> _num_samples_in_sample_set;
    std::vector<simd_t>     _subtree_sizes;

    void _compute(SetOfSampleSets const& samples) {
        KASSERT(_subtree_sizes.size() == 0ul, "Subtree sizes already computed.", tdt::assert::light);
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
            // Add all four bit counts bit-parallel
            _subtree_sizes[work_it->from()] += _subtree_sizes[work_it->to()];

            for (size_t sample_set_idx = 0; sample_set_idx < N; sample_set_idx++) {
                KASSERT(
                    (_subtree_sizes[work_it->from()][sample_set_idx]) <= _num_samples_in_sample_set[sample_set_idx],
                    "Number of samples below a node exceeds the number of samples in the tree sequence.",
                    tdt::assert::light
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
                    tdt::assert::light
                );
            }

            work_it++;
        }
    }
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

class NumSamplesBelowFactory {
public:
    template <typename BaseType = SampleId>
    static NumSamplesBelowAccessor<NumSamplesBelow<1, BaseType>>
    build(EdgeListGraph const& dag, SampleSet const& samples) {
        using SetOfSampleSets = typename NumSamplesBelow<1, BaseType>::SetOfSampleSets;
        auto const num_samples_below =
            std::make_shared<NumSamplesBelow<1, BaseType>>(dag, SetOfSampleSets{std::cref(samples)});
        return NumSamplesBelowAccessor(num_samples_below, 0);
    }

    template <typename BaseType = SampleId>
    static auto build(EdgeListGraph const& dag, SampleSet const& samples_0, SampleSet const& samples_1) {
        using SetOfSampleSets = typename NumSamplesBelow<2, BaseType>::SetOfSampleSets;
        auto num_samples_below =
            std::make_shared<NumSamplesBelow<2>>(dag, SetOfSampleSets{std::cref(samples_0), std::cref(samples_1)});
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
        using SetOfSampleSets  = typename NumSamplesBelow<4, BaseType>::SetOfSampleSets;
        auto num_samples_below = std::make_shared<NumSamplesBelow<4>>(
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
};
