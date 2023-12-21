#pragma once

#include <tskit/trees.h>

namespace sfkit::tskit {
class Children {
public:
    Children(tsk_tree_t const& tree, tsk_id_t const parent);
    class children_iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = tsk_id_t;
        using pointer           = value_type*;
        using reference         = value_type&;
        struct sentinel {};

        children_iterator(tsk_tree_t const& tree, tsk_id_t const parent);

        children_iterator& operator++();
        children_iterator  operator++(int);

        bool operator==(children_iterator const& other) const;
        bool operator!=(children_iterator const& other) const;
        bool operator==(sentinel);

        reference operator*();
        pointer   operator->();

    private:
        tsk_tree_t const& _tree;
        tsk_id_t          _child;
    };

    children_iterator           begin() const;
    children_iterator::sentinel end() const;

private:
    tsk_tree_t const& _tree;
    tsk_id_t const    _parent;
};
} // namespace sfkit::tskit
