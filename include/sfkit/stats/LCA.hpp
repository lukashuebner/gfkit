#pragma once

#include <array>
#include <vector>

#include <kassert/kassert.hpp>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/dag/DAGCompressedForest.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/samples/NumSamplesBelow.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/samples/primitives.hpp"
#include "sfkit/utils/checking_casts.hpp"

namespace sfkit::stats {

using sfkit::graph::EdgeListGraph;
using sfkit::graph::NodeId;
using sfkit::utils::asserting_cast;

class DAGLowestCommonAncestor {
public:
    DAGLowestCommonAncestor(EdgeListGraph const& dag) : _dag(dag) {
        // Check inputs
        KASSERT(_dag.check_postorder(), "DAG edges are not post-ordered.", sfkit::assert::normal);
        KASSERT(_dag.num_nodes() >= _dag.num_leaves(), "DAG has less nodes than leaves.", sfkit::assert::light);
    }

    [[nodiscard]] std::vector<NodeId> operator()(samples::SampleSet const& samples) const {
        return this->lca(samples);
    }

    [[nodiscard]] dag::SampleId num_nodes_in_dag() const {
        return asserting_cast<dag::SampleId>(_dag.num_nodes());
    }

    [[nodiscard]] dag::SampleId num_samples_in_dag() const {
        return asserting_cast<dag::SampleId>(_dag.num_leaves());
    }

    [[nodiscard]] std::vector<NodeId> lca(samples::SampleSet const& samples) const {
        std::vector<NodeId> lcas;

        KASSERT(
            _dag.num_leaves() <= samples.overall_num_samples(),
            "Number of leaves in the DAG is greater than he number of overall samples representable in the subtree "
            "size object. (NOT the number of samples actually in the SampleSet)",
            sfkit::assert::light
        );
        KASSERT(
            samples.popcount() <= samples.overall_num_samples(),
            "Number of samples in sample set is greater than he number of overall samples representable in the "
            "subtree size object. (NOT the number of samples actually in the SampleSet)",
            sfkit::assert::light
        );
        auto const num_requested_samples = samples.popcount();

        std::vector<lca_interm> subtree_sizes;
        subtree_sizes.resize(_dag.num_nodes(), {0, graph::INVALID_NODE_ID});
        for (dag::SampleId sample: samples) {
            subtree_sizes[sample] = {1, graph::INVALID_NODE_ID};
        }

        // Compute the subtree sizes using a post-order traversal
        auto edge_it = _dag.begin();
        while (edge_it != _dag.end()) {
            auto const& from = edge_it->from();
            auto const& to   = edge_it->to();

            if (subtree_sizes[from].lca != graph::INVALID_NODE_ID) {
                // LCA already propagated to this `from` node via another edge.
                KASSERT(subtree_sizes[from].samples_below == num_requested_samples);
                KASSERT(subtree_sizes[to].samples_below == 0u);
            } else {
                // LCA already below the `to` node
                if (subtree_sizes[to].samples_below == num_requested_samples) {
                    subtree_sizes[from].samples_below = num_requested_samples;
                    subtree_sizes[from].lca           = subtree_sizes[to].lca;
                } else { // LCA not below the `to` node
                    subtree_sizes[from].samples_below += subtree_sizes[to].samples_below;
                    // Is the current \c from node the LCA?
                    if (subtree_sizes[from].samples_below == num_requested_samples) {
                        subtree_sizes[from].lca = from;
                    }
                }
            }

            KASSERT(
                (subtree_sizes[edge_it->from()].samples_below) <= samples.overall_num_samples(),
                "Number of samples below a node exceeds the number of samples in the tree sequence.",
                sfkit::assert::light
            );

            edge_it++;
        }

        // Collect the per-tree LCAs
        lcas.reserve(_dag.num_trees());
        for (auto const root: _dag.roots()) {
            KASSERT(
                subtree_sizes[root].samples_below == samples.popcount(),
                "Number of samples below the root node does not match the number of samples in the sample set.",
                sfkit::assert::light
            );
            lcas.push_back(subtree_sizes[root].lca);
        }
        return lcas;
    }

private:
    EdgeListGraph const& _dag; // As a post-order sorted edge list

    struct lca_interm {
        NodeId samples_below;
        NodeId lca;
    };
};

} // namespace sfkit::stats
