#include <sfkit/tskit/ChangedNodesView.hpp>
#include <tskit/trees.h>

namespace sfkit::tskit {

class TSKitTreeSequence;

ChangedNodesView::ChangedNodesView(
    tsk_tree_t const& tree, TSKitTreeSequence const& tree_sequence, tsk_tree_position_t const& tree_pos
)
    : _tree(tree),
      _tree_sequence(tree_sequence),
      _tree_pos(tree_pos) {}

} // namespace sfkit::tskit
