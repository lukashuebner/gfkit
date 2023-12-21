#include "sfkit/tskit/EulertourView.hpp"

#include <kassert/kassert.hpp>
#include <tskit/trees.h>

#include "sfkit/tskit/TSKitTreeSequence.hpp"

namespace sfkit::tskit {

EulertourView::EulertourView(tsk_tree_t const& tree, TSKitTreeSequence const& tree_sequence)
    : _tree(tree),
      _tree_sequence(tree_sequence) {}

EulertourView::iterator::iterator(tsk_tree_t const& tree, TSKitTreeSequence const& tree_sequence)
    : _tree(tree),
      _just_moved_up(false),
      _tree_sequence(tree_sequence) {
    KASSERT(_tree.virtual_root != TSK_NULL, "The virtual root of the tree does not exist.", sfkit::assert::light);
    KASSERT(
        _tree.right_child[_tree.virtual_root] != TSK_NULL,
        "The virtual root of the tree has no children.",
        sfkit::assert::light
    );
    KASSERT(
        _tree.left_child[_tree.virtual_root] == _tree.right_child[_tree.virtual_root],
        "The virtual root of the tree has more than one child.",
        sfkit::assert::light
    );
    _current_node = _tree.left_child[_tree.virtual_root];
}

EulertourView::iterator& EulertourView::iterator::operator++() {
    if (_just_moved_up || is_sample()) {
        // The children of this node were already processed OR there are no children because the current
        // node is a sample node.
        if (_tree.right_sib[_current_node] != TSK_NULL) {
            // Move to the next unprocessed sibling of this node.
            _current_node  = _tree.right_sib[_current_node];
            _just_moved_up = false;
        } else {
            // There are no more siblings of this node -> move up.
            _current_node  = _tree.parent[_current_node];
            _just_moved_up = true;
        }
    } else {
        // We just moved down or right (to a sibling), the children of this node are not processed yet
        //   -> move down
        _current_node = _tree.left_child[_current_node];
    }
    return *this;
}

EulertourView::iterator EulertourView::iterator::operator++(int) {
    iterator tmp(*this);

    this->operator++();
    return tmp;
}

bool EulertourView::iterator::operator==(iterator const& other) const {
    return _current_node == other._current_node && &_tree == &other._tree && _just_moved_up == other._just_moved_up;
}

bool EulertourView::iterator::operator!=(iterator const& other) const {
    return !(*this == other);
}

bool EulertourView::iterator::operator==(sentinel) {
    return _current_node == TSK_NULL;
}

EulertourView::iterator::reference EulertourView::iterator::operator*() {
    KASSERT(_current_node != TSK_NULL);
    return _current_node;
}

EulertourView::iterator::reference EulertourView::iterator::node_id() {
    return operator*();
}

bool EulertourView::iterator::first_visit() {
    return !_just_moved_up;
}

bool EulertourView::iterator::second_visit() {
    return _just_moved_up;
}

bool EulertourView::iterator::is_sample() {
    KASSERT(_current_node != TSK_NULL);
    return _tree_sequence.is_sample(_current_node);
}

EulertourView::iterator::pointer EulertourView::iterator::operator->() {
    KASSERT(_current_node != TSK_NULL);
    return &_current_node;
}

EulertourView::iterator EulertourView::begin() const {
    return iterator{_tree, _tree_sequence};
}

EulertourView::iterator::sentinel EulertourView::end() const {
    return iterator::sentinel{};
}

} // namespace sfkit::tskit
