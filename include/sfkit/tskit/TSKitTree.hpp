#pragma once

#include <cstddef>
#include <cstdint>
#include <ranges>
#include <unordered_set>
#include <vector>

#include <tskit/trees.h>

#include "sfkit/graph/primitives.hpp"
#include "sfkit/include-redirects/hopscotch_set.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/samples/primitives.hpp"
#include "sfkit/sequence/Mutation.hpp"
#include "sfkit/sequence/Sequence.hpp"
#include "sfkit/tskit/ChangedNodesView.hpp"
#include "sfkit/tskit/ChildrenView.hpp"
#include "sfkit/tskit/EulertourView.hpp"

namespace sfkit::tskit {

constexpr int TSK_NULL_TREE = 0;

using sfkit::graph::EdgeId;
using sfkit::graph::NodeId;
using sfkit::graph::TreeId;
using sfkit::samples::SampleId;
using sfkit::sequence::MutationId;
using sfkit::sequence::SiteId;

class TSKitTreeSequence;

// TODO Change to sfkit:tskit::Tree
class TSKitTree {
public:
    TSKitTree(TSKitTreeSequence& tree_sequence);
    ~TSKitTree();
    bool first();
    bool next();

    [[nodiscard]] bool is_valid() const;
    [[nodiscard]] bool is_tree() const;
    [[nodiscard]] bool is_null() const;

    [[nodiscard]] tsk_id_t    tree_id() const;
    [[nodiscard]] std::size_t num_roots() const;
    [[nodiscard]] std::size_t num_samples() const;
    [[nodiscard]] tsk_id_t    num_children(tsk_id_t node) const;
    [[nodiscard]] std::size_t max_node_id() const;
    [[nodiscard]] tsk_id_t    root() const;
    [[nodiscard]] bool        is_root(tsk_id_t const node) const;
    [[nodiscard]] bool        is_sample(tsk_id_t node) const;
    [[nodiscard]] NodeId      lca(tsk_id_t const u, tsk_id_t const v);

    [[nodiscard]] std::span<tsk_id_t> postorder();
    [[nodiscard]] EulertourView       eulertour() const;
    // [[nodiscard]] ChangedNodesView    changed_nodes();
    [[nodiscard]] Children children(tsk_id_t const parent);
    [[nodiscard]] tsk_id_t parent(tsk_id_t const node) const;

    tsl::hopscotch_set<tsk_id_t> invalidated_nodes() const;

private:
    [[nodiscard]] size_t _current_tree_size_bound() const;
    void                 _postorder_nodes_resize();
    // TODO Create a variadic templated _create_tsk_sample_sets and use it in the functions above.

    TSKitTreeSequence&    _tree_sequence;
    tsk_tree_t            _tree;
    tsk_tree_position_t   _tree_pos;
    tsk_id_t              _tree_id;
    int                   _state = -1;
    std::vector<tsk_id_t> _postorder_nodes;
};

} // namespace sfkit::tskit
