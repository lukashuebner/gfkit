#include "sfkit/tskit/ChildrenView.hpp"

namespace sfkit::tskit {

// TODO Rename to ChildrenView
Children::Children(tsk_tree_t const& tree, tsk_id_t const parent) : _tree(tree), _parent(parent) {}

Children::children_iterator::children_iterator(tsk_tree_t const& tree, tsk_id_t const parent)
    : _tree(tree),
      _child(tree.left_child[parent]) {}

Children::children_iterator& Children::children_iterator::operator++() {
    _child = _tree.right_sib[_child];
    // This provided a small but measurable speedup (~10%)
    __builtin_prefetch(&_tree.right_sib[_child]);
    return *this;
}

Children::children_iterator Children::children_iterator::operator++(int) {
    children_iterator tmp(*this);
    this->            operator++();
    return tmp;
}

bool Children::children_iterator::operator==(children_iterator const& other) const {
    return _child == other._child && &_tree == &other._tree;
}

bool Children::children_iterator::operator!=(children_iterator const& other) const {
    return !(*this == other);
}

bool Children::children_iterator::operator==(sentinel) {
    return _child == TSK_NULL;
}

Children::children_iterator::reference Children::children_iterator::operator*() {
    return _child;
}

Children::children_iterator::pointer Children::children_iterator::operator->() {
    return &_child;
}

Children::children_iterator Children::begin() const {
    return children_iterator{_tree, _parent};
}

Children::children_iterator::sentinel Children::end() const {
    return children_iterator::sentinel{};
}

} // namespace sfkit::tskit
