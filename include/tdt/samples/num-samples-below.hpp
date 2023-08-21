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

    [[nodiscard]] SampleId num_samples_below(NodeId node_id) const {
        KASSERT(node_id < _subtree_sizes.size(), "Subtree ID out of bounds.", tdt::assert::light);
        return _subtree_sizes[node_id];
    }

    [[nodiscard]] SampleId operator()(NodeId node_id) const {
        return this->num_samples_below(node_id);
    }

    [[nodiscard]] SampleId operator[](NodeId node_id) const {
        return this->num_samples_below(node_id);
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
            _dag.num_leaves() <= samples.overall_num_samples(),
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

        constexpr auto num_edges_to_prefetch = 128;
        auto           prefetch_it           = _dag.begin();
        for (size_t i = 0; i < num_edges_to_prefetch && prefetch_it != _dag.end(); i++) {
            __builtin_prefetch(&_subtree_sizes[prefetch_it->from()], 0, 3);
            __builtin_prefetch(&_subtree_sizes[prefetch_it->to()], 0, 3);
            prefetch_it++;
        }

        auto work_it = _dag.begin();
        while (prefetch_it != _dag.end()) {
            _subtree_sizes[work_it->from()] += _subtree_sizes[work_it->to()];
            KASSERT(
                _subtree_sizes[work_it->from()] <= samples.popcount(),
                "Number of samples below a node exceeds the number of samples in the tree sequence.",
                tdt::assert::light
            );
            work_it++;
            __builtin_prefetch(&_subtree_sizes[prefetch_it->from()], 0, 3);
            __builtin_prefetch(&_subtree_sizes[prefetch_it->to()], 0, 3);
            prefetch_it++;
        }

        while (work_it != _dag.end()) {
            _subtree_sizes[work_it->from()] += _subtree_sizes[work_it->to()];
            KASSERT(
                _subtree_sizes[work_it->from()] <= samples.popcount(),
                "Number of samples below a node exceeds the number of samples in the tree sequence.",
                tdt::assert::light
            );
            work_it++;
        }
    }
};
