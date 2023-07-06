#pragma once

#include <vector>

#include <kassert/kassert.hpp>

#include "tdt/graph/edge-list-graph.hpp"
// TODO Unify naming of files
#include "tdt/assertion_levels.hpp"
#include "tdt/samples/sample-set.hpp"

class NumSamplesBelow {
public:
    NumSamplesBelow(EdgeListGraph const& dag_postorder_edges, SampleSet const& samples) : _dag(dag_postorder_edges) {
        compute(samples);
    }

    SampleId num_samples_below(SubtreeId subtree_id) {
        KASSERT(subtree_id < _subtree_sizes.size(), "Subtree ID out of bounds.", tdt::assert::light);
        return _subtree_sizes[subtree_id];
    }

    SampleId operator[](SubtreeId subtree_id) {
        return this->num_samples_below(subtree_id);
    }

    SampleId num_nodes() const {
        return asserting_cast<SampleId>(_dag.num_nodes());
    }

private:
    EdgeListGraph const&  _dag; // As a post-order sorted edge list
    std::vector<SampleId> _subtree_sizes;

    void compute(SampleSet const& samples) {
        KASSERT(_dag.check_postorder(), "DAG edges are not post-ordered.", tdt::assert::normal);

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
