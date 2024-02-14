#pragma once

#include <tskit/trees.h>

#include "sfkit/tskit/TSKitTreeSequence.hpp"

namespace sfkit::tskit {

class TSKitTreeSequence;

// --- ChangedNodesView ---
// class ChangedNodesView {
// public:
//     ChangedNodesView(
//         tsk_tree_t const& tree, TSKitTreeSequence const& tree_sequence, tsk_tree_position_t const& tree_pos
//     );

//     class iterator {
//     public:
//         using iterator_category = std::forward_iterator_tag;
//         using difference_type   = std::ptrdiff_t;
//         using value_type        = tsk_id_t;
//         using pointer           = value_type*;
//         using reference         = value_type&;
//         struct sentinel {};

//         iterator(){};
//     };

// private:
//     tsk_tree_t const&          _tree;
//     TSKitTreeSequence const&   _tree_sequence;
//     tsk_tree_position_t const& _tree_pos;
// };

} // namespace sfkit::tskit
