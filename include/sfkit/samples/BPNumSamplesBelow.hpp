#pragma once

#include "plf_stack.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <experimental/simd>

#include "sfkit/bp/BPCompressedForest.hpp"
#include "sfkit/samples/BPNumSamplesBelow.hpp"
#include "sfkit/samples/DAGNumSamplesBelow.hpp"
#include "sfkit/samples/NumSamplesBelow.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/samples/primitives.hpp"
#include "sfkit/utils/BufferedSDSLBitVectorView.hpp"

namespace sfkit::samples {

using sfkit::bp::BPCompressedForest;
using sfkit::graph::NodeId;
using sfkit::samples::SampleId;
using sfkit::utils::asserting_cast;
using sfkit::utils::BufferedSDSLBitVectorView;
namespace stdx = std::experimental;

template <size_t N, typename BaseType>
class NumSamplesBelow<BPCompressedForest, N, BaseType> {
public:
    NumSamplesBelow(BPCompressedForest const& forest, SetOfSampleSets<N> const& samples) : _forest(forest) {
        // Check inputs
        KASSERT(_forest.num_nodes() >= _forest.num_leaves(), "DAG has less nodes than leaves.", sfkit::assert::light);
        for (auto const& sample_set: samples) {
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
            for (auto const& other_sample_set: samples) {
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

    void _compute(SetOfSampleSets<N> const& samples) {
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

        std::vector<simd_t>  sample_counts;
        plf::stack<SampleId> num_children;
        sample_counts.resize(_forest.num_samples() + 1); // TODO Think about the maximum size
        // We're wasting the first entry here as a sentinel to simplify the logic below.
        auto sample_counts_top = sample_counts.begin();

        NodeId   inner_node_id = _forest.num_samples();
        SampleId leaf_rank     = 0;
        size_t   ref_rank      = 0;

        // TODO If the buffered iterator is not faster, remove it again
        KASSERT(bp.size() == is_ref.size(), "balanced_parenthesis and is_reference are of different size");
        BufferedSDSLBitVectorView bp_view{bp};
        BufferedSDSLBitVectorView is_ref_view{is_ref};

        auto bp_it      = bp_view.begin();
        auto bp_end     = bp_view.end();
        auto is_ref_it  = is_ref_view.begin();
        auto is_ref_end = is_ref_view.end();

        bool   last_bp = bp::PARENS_CLOSE;
        size_t level   = 0; // Distance from root
        // bool last_is_leaf = false;
        num_children.emplace(0);
        while (bp_it != bp_end) {
            KASSERT(num_children.size() == level + 1);
            KASSERT(bp_it != bp_end && is_ref_it != is_ref_end);
            if (*is_ref_it) { // reference
                // Reference always occur as tuples of (open, close) in BP and (true, true) in is_ref
                KASSERT(*bp_it == bp::PARENS_OPEN);
                ++bp_it;
                ++is_ref_it;
                KASSERT(bp_it != bp_end && is_ref_it != is_ref_end);
                KASSERT(*is_ref_it);
                KASSERT(*bp_it == bp::PARENS_CLOSE);

                NodeId const node_id = _forest.node_id_ref_by_rank(ref_rank);
                if (level > 0) [[likely]] { // We're not referring to a whole tree
                    ++sample_counts_top;
                    *sample_counts_top = _subtree_sizes[node_id];
                    num_children.top()++;
                }
                ++ref_rank;
            } else { // description of subtree
                if (*bp_it == bp::PARENS_OPEN) {
                    ++level;
                    num_children.top()++;
                    num_children.emplace(0);
                } else { // &bp_it == bp::PARENS_CLOSE
                    --level;
                    if (last_bp == bp::PARENS_OPEN) { // sample
                        KASSERT(leaf_rank < _forest.num_samples());
                        SampleId const leaf_id = _forest.leaf_idx_to_id(leaf_rank);
                        ++sample_counts_top;
                        *sample_counts_top = _subtree_sizes[leaf_id];
                        ++leaf_rank;
                    } else { // inner node
                        KASSERT(sample_counts.size() >= 2ul);
                        for ([[maybe_unused]] SampleId child = 0; child < num_children.top() - 1; ++child) {
                            auto const other = *sample_counts_top;
                            --sample_counts_top;
                            *sample_counts_top += other;
                        }
                        KASSERT(inner_node_id < _forest.num_nodes());
                        _subtree_sizes[inner_node_id] = *sample_counts_top;
                        if (level == 0) [[unlikely]] {
                            --sample_counts_top;
                        }
                        ++inner_node_id;
                    }
                    num_children.pop();
                }
            }
            last_bp = *bp_it;
            ++bp_it;
            ++is_ref_it;
        }
        KASSERT(sample_counts_top - sample_counts.begin() == 0);
        KASSERT(num_children.size() == 1ul);
    }
};
} // namespace sfkit::samples
