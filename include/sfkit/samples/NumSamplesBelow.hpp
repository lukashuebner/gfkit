#pragma once

#include <array>
#include <experimental/simd>
#include <vector>

#include <kassert/kassert.hpp>
#include <plf_stack.h>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/graph/BPCompressedForest.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/utils/BufferedSDSLBitVectorView.hpp"

namespace stdx = std::experimental;
using sfkit::utils::BufferedSDSLBitVectorView;

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
            // Add all four bit counts bit-parallel
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

// TODO Split up this and other files.
// TODO Move from header-only to .hpp + .cpp
template <size_t N = 1, typename BaseType = SampleId>
class BPNumSamplesBelow {
public:
    using SetOfSampleSets = std::array<std::reference_wrapper<SampleSet const>, N>;

    BPNumSamplesBelow(BPCompressedForest const& forest, SetOfSampleSets const& samples) : _forest(forest) {
        // Check inputs
        KASSERT(_forest.num_nodes() >= _forest.num_leaves(), "DAG has less nodes than leaves.", sfkit::assert::light);
        for (auto sample_set: samples) {
            KASSERT(
                _forest.num_leaves() <= sample_set.get().overall_num_samples(),
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
                    "actually in the SampleSet)",
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
        return asserting_cast<SampleId>(_forest.num_nodes());
    }

    [[nodiscard]] SampleId num_samples_in_dag() const {
        return asserting_cast<SampleId>(_forest.num_leaves());
    }

    [[nodiscard]] SampleId num_samples_in_sample_set(SampleSetId sample_set_id) const {
        KASSERT(sample_set_id >= 0 && sample_set_id <= N, "Sample set ID invalid.", sfkit::assert::light);
        return _num_samples_in_sample_set[sample_set_id];
    }

private:
    using simd_t = stdx::fixed_size_simd<BaseType, N>;
    BPCompressedForest const& _forest;
    std::array<SampleId, N>   _num_samples_in_sample_set;
    std::vector<simd_t>       _subtree_sizes;

    void _compute(SetOfSampleSets const& samples) {
        KASSERT(_subtree_sizes.size() == 0ul, "Subtree sizes already computed.", sfkit::assert::light);
        _subtree_sizes.resize(_forest.num_nodes(), 0);

        KASSERT(samples.size() == N);
        for (size_t sample_set_idx = 0; sample_set_idx < N; sample_set_idx++) {
            for (SampleId sample: samples[sample_set_idx].get()) {
                _subtree_sizes[sample][sample_set_idx] = 1;
            }
        }

        auto const& bp     = _forest.balanced_parenthesis();
        auto const& is_ref = _forest.is_reference();
        KASSERT(
            bp.size() == is_ref.size(),
            "balanced_parenthesis and is_reference are of different size",
            sfkit::assert::light
        );

        KASSERT(samples.size() == N);
        for (size_t sample_set_idx = 0; sample_set_idx < N; sample_set_idx++) {
            for (SampleId sample: samples[sample_set_idx].get()) {
                _subtree_sizes[sample][sample_set_idx] = 1;
            }
        }

        // TODO Abstract this away into a NodeIdCounter class which we also use during forest compression
        plf::stack<simd_t> sample_counts;
        sample_counts.reserve(_forest.num_samples());

        NodeId   inner_node_id = _forest.num_samples();
        SampleId leaf_rank     = 0;
        size_t   ref_rank      = 0;

        // TODO It the buffered iterator is not faster, remove it again
        // auto bp_size = bp.size(); // Doing this in the for-loop used a substantial amount of runtime.
        KASSERT(bp.size() == is_ref.size(), "balanced_parenthesis and is_reference are of different size");
        // TODO Write a faster iterator for the BP sequence (extracting 64bit, shifting over them...)
        auto bp_view     = BufferedSDSLBitVectorView(bp);
        auto bp_it       = bp_view.begin();
        auto bp_end      = bp_view.end();
        auto is_ref_view = BufferedSDSLBitVectorView(is_ref);
        auto is_ref_it   = is_ref_view.begin();
        auto is_ref_end  = is_ref_view.end();
        // auto is_leaf_view = BufferedSDSLBitVectorView(_forest.is_leaf());
        // auto is_leaf_it = is_leaf_view.begin();
        // auto is_leaf_end = is_leaf_view.end();
        // TODO Idea: use less bits to represent is_ref, is_leaf, bp and interleave the arrays for more efficient memory
        // access.

        bool last_bp = bp::PARENS_CLOSE;
        // bool last_is_leaf = false;
        while (bp_it != bp_end) {
            // TODO Add operator!= for sentinel
            KASSERT(!(bp_it == bp_end) && !(is_ref_it == is_ref_end));
            if (*is_ref_it) { // reference
                // Reference always occur as tuples of (open, close) in BP and (true, true) in is_ref
                KASSERT(*bp_it == bp::PARENS_OPEN);
                ++bp_it;
                ++is_ref_it;
                // TODO Change back
                KASSERT(!(bp_it == bp_end) && !(is_ref_it == is_ref_end));
                // KASSERT(bp_it != bp_end && is_ref_it != is_ref_end);
                KASSERT(*is_ref_it);
                KASSERT(*bp_it == bp::PARENS_CLOSE);

                NodeId const node_id = _forest.node_id_ref_by_rank(ref_rank);
                sample_counts.push(_subtree_sizes[node_id]);
                ++ref_rank;
            } else { // description of subtree
                if (*bp_it == bp::PARENS_CLOSE) {
                    // KASSERT(idx >= 1ul); // The first bit in the sequence can never be a closing parenthesis.
                    if (last_bp == bp::PARENS_OPEN) { // sample
                        // KASSERT(*is_leaf_it);
                        // KASSERT(last_is_leaf);
                        KASSERT(leaf_rank < _forest.num_samples());
                        SampleId const leaf_id = _forest.leaf_idx_to_id(leaf_rank);
                        // KASSERT(leaf_id == _forest.node_id(idx - 1));
                        sample_counts.push(_subtree_sizes[leaf_id]);
                        ++leaf_rank;
                    } else { // inner node
                        KASSERT(sample_counts.size() >= 2ul);
                        auto const other = sample_counts.top();
                        sample_counts.pop();
                        sample_counts.top() += other;
                        KASSERT(inner_node_id < _forest.num_nodes());
                        // KASSERT(_forest.node_id(asserting_cast<size_t>(idx)) == inner_node_id);
                        _subtree_sizes[inner_node_id] = sample_counts.top();
                        ++inner_node_id;
                    }
                }
            }
            last_bp = *bp_it;
            // last_is_leaf = *is_leaf_it;
            ++bp_it;
            ++is_ref_it;
            // is_leaf_it++;
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
        using SetOfSampleSets          = typename NumSamplesBelow<1, BaseType>::SetOfSampleSets;
        constexpr auto num_sample_sets = 1;

        auto const num_samples_below =
            std::make_shared<NumSamplesBelow<num_sample_sets, BaseType>>(dag, SetOfSampleSets{std::cref(samples)});

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
        using SetOfSampleSets          = typename NumSamplesBelow<2, BaseType>::SetOfSampleSets;
        constexpr auto num_sample_sets = 2;

        auto num_samples_below = std::make_shared<NumSamplesBelow<num_sample_sets>>(
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
        using SetOfSampleSets          = typename NumSamplesBelow<4, BaseType>::SetOfSampleSets;
        constexpr auto num_sample_sets = 4;

        auto num_samples_below = std::make_shared<NumSamplesBelow<num_sample_sets, BaseType>>(
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
