
#include "sfkit/tskit/TSKitTree.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <unordered_set>
#include <vector>

#include <kassert/kassert.hpp>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/include-redirects/hopscotch_set.hpp"
#include "sfkit/tskit/ChangedNodesView.hpp"
#include "sfkit/tskit/ChildrenView.hpp"
#include "sfkit/tskit/EulertourView.hpp"
#include "sfkit/tskit/MutationsView.hpp"
#include "sfkit/tskit/TSKitTreeSequence.hpp"

namespace sfkit::tskit {

using sfkit::utils::asserting_cast;

TSKitTree::TSKitTree(TSKitTreeSequence& tree_sequence) : _tree_sequence(tree_sequence) {
    int ret = tsk_tree_init(&_tree, &(_tree_sequence.underlying()), 0);
    KASSERT(ret == 0, "Failed to initialize the tree data structure", sfkit::assert::light);

    ret = tsk_tree_position_init(&_tree_pos, &_tree_sequence.underlying(), 0);
    KASSERT(ret == 0, "Failed to initialize the tree position.", sfkit::assert::light);

    // Load the first tree
    first();
}

TSKitTree::~TSKitTree() {
    tsk_tree_free(&_tree);
    tsk_tree_position_free(&_tree_pos);
}

bool TSKitTree::first() {
    _tree_id = 0;

    _state = tsk_tree_first(&_tree);
    KASSERT(is_valid(), "Failed to goto the first tree.", sfkit::assert::light);

    // Is there no other way to reset the tree position?
    // TODO Use std::optional to avoid double allocation when creating.
    tsk_tree_position_free(&_tree_pos);
    int const ret = tsk_tree_position_init(&_tree_pos, &_tree_sequence.underlying(), 0);
    KASSERT(ret == 0, "Failed to reset the tree position.", sfkit::assert::light);
    if (is_tree()) {
        tsk_tree_position_next(&_tree_pos);
    }

    return is_tree();
}

bool TSKitTree::next() {
    ++_tree_id;

    _state = tsk_tree_next(&_tree);
    KASSERT(is_valid(), "Failed to goto the next tree.", sfkit::assert::light);

    tsk_tree_position_next(&_tree_pos);

    return is_tree();
}

tsk_id_t TSKitTree::tree_id() const {
    return _tree_id;
}

bool TSKitTree::is_valid() const {
    return is_tree() || is_null();
}

bool TSKitTree::is_tree() const {
    return _state == TSK_TREE_OK;
}

bool TSKitTree::is_null() const {
    return _state == TSK_NULL_TREE;
}

std::size_t TSKitTree::num_roots() const {
    KASSERT(is_valid(), "The tree is not valid.", sfkit::assert::light);
    return tsk_tree_get_num_roots(&_tree);
}

std::size_t TSKitTree::max_node_id() const {
    return _tree_sequence.num_nodes();
}

std::size_t TSKitTree::num_samples() const {
    return _tree_sequence.num_samples();
}

bool TSKitTree::is_sample(tsk_id_t node) const {
    KASSERT(is_valid(), "The tree is not valid.", sfkit::assert::light);
    KASSERT(asserting_cast<tsk_size_t>(node) <= max_node_id(), "The node is not valid.", sfkit::assert::light);
    return _tree_sequence.is_sample(node);
}

tsk_id_t TSKitTree::num_children(tsk_id_t node) const {
    KASSERT(is_valid(), "The tree is not valid.", sfkit::assert::light);
    KASSERT(asserting_cast<tsk_size_t>(node) <= max_node_id(), "The node is not valid.", sfkit::assert::light);
    return _tree.num_children[node];
}

tsk_id_t TSKitTree::root() const {
    KASSERT(is_valid(), "The tree is not valid.", sfkit::assert::light);
    KASSERT(
        tsk_tree_get_num_roots(&_tree) == 1ul,
        "The tskit tree has more than one root. We do not support this yet."
    );
    KASSERT(_tree.virtual_root != TSK_NULL);
    KASSERT(_tree.right_child[_tree.virtual_root] != TSK_NULL);
    KASSERT(_tree.left_child[_tree.virtual_root] == _tree.right_child[_tree.virtual_root]);
    KASSERT(num_children(_tree.virtual_root) == 1);
    return _tree.right_child[_tree.virtual_root];
}

bool TSKitTree::is_root(tsk_id_t const node) const {
    // The following would NOT work, as if the node is not in the current tree, the parent is also set to TKS_NULL!
    // return _tree.parent[node] == TSK_NULL;

    return root() == node;
}

std::span<tsk_id_t> TSKitTree::postorder() {
    int        ret;
    tsk_size_t num_nodes;
    _postorder_nodes_resize();
    ret = tsk_tree_postorder(&_tree, _postorder_nodes.data(), &num_nodes);
    KASSERT(ret == 0, "Failed to get the postorder traversal of the tree.", sfkit::assert::light);
    KASSERT(
        num_nodes <= _postorder_nodes.size(),
        "The number of nodes in the postorder traversal is too large to fit into the buffer.",
        sfkit::assert::light
    );
    return std::span{_postorder_nodes}.subspan(0, num_nodes);
}

tsl::hopscotch_set<tsk_id_t> TSKitTree::invalidated_nodes() const {
    tsl::hopscotch_set<tsk_id_t> changed_nodes;

    tsk_id_t const*       in_it   = _tree_pos.in.order + _tree_pos.in.start;
    tsk_id_t const* const in_end  = _tree_pos.in.order + _tree_pos.in.stop;
    tsk_id_t const*       out_it  = _tree_pos.out.order + _tree_pos.out.start;
    tsk_id_t const* const out_end = _tree_pos.out.order + _tree_pos.out.stop;

    for (auto idx = _tree_pos.out.start; idx < _tree_pos.out.stop; ++idx) {
        auto const e      = _tree_pos.out.order[idx];
        auto const parent = _tree_pos.tree_sequence->tables->edges.parent[e];
        changed_nodes.insert(parent);
    }

    for (auto idx = _tree_pos.in.start; idx < _tree_pos.in.stop; ++idx) {
        auto const e      = _tree_pos.in.order[idx];
        auto const parent = _tree_pos.tree_sequence->tables->edges.parent[e];
        changed_nodes.insert(parent);
    }

    // TODO Assert there are no leaves

    // Propagate the 'changed' tag upwards the tree
    std::vector<tsk_id_t> queue;
    queue.reserve(_tree.num_nodes);
    std::copy(changed_nodes.begin(), changed_nodes.end(), std::back_inserter(queue));

    auto node_it = queue.begin();
    while (node_it != queue.end()) {
        auto const node   = *node_it;
        auto const parent = _tree.parent[node];

        if (parent != TSK_NULL && !changed_nodes.contains(parent)) {
            KASSERT(
                queue.capacity() > queue.size() + 1,
                "The queue would need to be resized, resulting in iterator invalidation.",
                sfkit::assert::light
            );
            queue.push_back(parent);
            changed_nodes.insert(parent);
        }

        ++node_it;
    }

    return changed_nodes;
}

EulertourView TSKitTree::eulertour() const {
    return EulertourView{_tree, _tree_sequence};
}

Children TSKitTree::children(tsk_id_t const parent) {
    return Children{_tree, parent};
}

tsk_id_t TSKitTree::parent(tsk_id_t const node) const {
    KASSERT(is_tree(), "The tree is not valid.", sfkit::assert::light);
    KASSERT(asserting_cast<tsk_size_t>(node) <= max_node_id(), "The node is not valid.", sfkit::assert::light);
    return _tree.parent[node];
}

size_t TSKitTree::_current_tree_size_bound() const {
    return tsk_tree_get_size_bound(&_tree);
}

// Reserve space for the traversal lists. Sadly, tsk provides us only an upper bound /on this tree/.
// We therefore might need to resize the vectors later on.
void TSKitTree::_postorder_nodes_resize() {
    if (_postorder_nodes.size() < _current_tree_size_bound()) { // We never want to shrink the vector
        _postorder_nodes.resize(_current_tree_size_bound(), TSK_NULL);
    }
}

} // namespace sfkit::tskit
