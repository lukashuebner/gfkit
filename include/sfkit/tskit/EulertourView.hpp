#pragma once

#include <tskit/trees.h>

namespace sfkit::tskit {

class TSKitTreeSequence;

class EulertourView {
public:
    EulertourView(tsk_tree_t const& tree, TSKitTreeSequence const& tree_sequence);

    class iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = tsk_id_t;
        using pointer           = value_type*;
        using reference         = value_type&;
        struct sentinel {};

        iterator(tsk_tree_t const& tree, TSKitTreeSequence const& tree_sequence);
        iterator&               operator++();
        iterator                operator++(int);
        [[nodiscard]] bool      operator==(iterator const& other) const;
        [[nodiscard]] bool      operator!=(iterator const& other) const;
        [[nodiscard]] bool      operator==(sentinel);
        [[nodiscard]] reference operator*();
        [[nodiscard]] reference node_id();
        [[nodiscard]] bool      first_visit();
        [[nodiscard]] bool      second_visit();
        [[nodiscard]] bool      is_sample();
        [[nodiscard]] pointer   operator->();

    private:
        tsk_tree_t const&        _tree;
        tsk_id_t                 _current_node;
        bool                     _just_moved_up;
        TSKitTreeSequence const& _tree_sequence;
    };

    iterator           begin() const;
    iterator::sentinel end() const;

private:
    tsk_tree_t const&        _tree;
    TSKitTreeSequence const& _tree_sequence;
};

} // namespace sfkit::tskit
