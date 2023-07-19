#pragma once

#include <vector>

#include <kassert/kassert.hpp>

#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/samples/sample-set.hpp"
// TODO Unify naming of files
#include "tdt/assertion_levels.hpp"
#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/samples/sample-set.hpp"

class NumSamplesBelow {
public:
    NumSamplesBelow(EdgeListGraph const& dag_postorder_edges, SampleSet const& samples)
        : _dag(dag_postorder_edges),
          _num_samples_in_sample_set(samples.popcount()) {
        compute(samples);
    }

    [[nodiscard]] SampleId num_samples_below(SubtreeId subtree_id) const {
        KASSERT(subtree_id < _subtree_sizes.size(), "Subtree ID out of bounds.", tdt::assert::light);
        return _subtree_sizes[subtree_id];
    }

    [[nodiscard]] SampleId operator()(SubtreeId subtree_id) const {
        return this->num_samples_below(subtree_id);
    }

    [[nodiscard]] SampleId operator[](SubtreeId subtree_id) const {
        return this->num_samples_below(subtree_id);
    }

    [[nodiscard]] SampleId num_nodes_in_dag() const {
        return asserting_cast<SampleId>(_dag.num_nodes());
    }

    [[nodiscard]] SampleId num_samples_in_dag() const {
        return asserting_cast<SampleId>(_dag.num_leaves());
    }

    [[nodiscard]] SampleId num_samples_in_sample_set() const {
        return _num_samples_in_sample_set;
    }

private:
    EdgeListGraph const&  _dag; // As a post-order sorted edge list
    SampleId              _num_samples_in_sample_set;
    std::vector<SampleId> _subtree_sizes;

    void compute(SampleSet const& samples) {
        KASSERT(_dag.check_postorder(), "DAG edges are not post-ordered.", tdt::assert::normal);
        KASSERT(
            _dag.num_leaves() <= samples.num_nodes_in_dag(),
            "Number of leaves in the DAG is greater than he number of overall samples representable in the subtree "
            "size object. (NOT the number of samples actually in the SampleSet)",
            tdt::assert::light
        );

        KASSERT(_dag.num_nodes() >= _dag.num_leaves(), "DAG has less nodes than leaves.", tdt::assert::light);
        KASSERT(_subtree_sizes.size() == 0ul, "Subtree sizes already computed.", tdt::assert::light);
        _subtree_sizes.resize(_dag.num_nodes(), 0);
        for (auto&& sample: samples) {
            _subtree_sizes[sample] = 1;
        }

        // TODO Try prefetching a few edges ahead
        for (auto&& edge: _dag) {
            _subtree_sizes[edge.from()] += _subtree_sizes[edge.to()];
            KASSERT(
                _subtree_sizes[edge.from()] <= samples.popcount(),
                "Number of samples below a node exceeds the number of samples in the tree sequence.",
                tdt::assert::light
            );
        }
    }
};
