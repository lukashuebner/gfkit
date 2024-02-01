#pragma once

#include <cstddef>
#include <iterator>
#include <numeric>
#include <span>

#include <kassert/kassert.hpp>
#include <tskit.h>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/graph/primitives.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/sequence/Mutation.hpp"
#include "sfkit/sequence/Sequence.hpp"
#include "sfkit/utils/checking_casts.hpp"

// namespace sfkit::tskit {

// constexpr int TSK_NULL_TREE = 0;

// using sfkit::graph::EdgeId;
// using sfkit::graph::NodeId;
// using sfkit::graph::TreeId;
// using sfkit::samples::SampleId;
// using sfkit::samples::SampleSet;
// using sfkit::sequence::MutationId;
// using sfkit::sequence::SiteId;
// using sfkit::utils::asserting_cast;

// using TskMutationView = std::span<tsk_mutation_t const>;

// class TSKitTreeSequence {
// public:
//     TSKitTreeSequence(std::string const& trees_file) : _trees_file(trees_file) {
//         int ret;

//         // Load the tree sequence from the .trees file
//         ret = tsk_treeseq_load(&_tree_sequence, _trees_file.c_str(), 0);
//         KASSERT(tskit_noerr(ret), "Failed to load tree sequence from the .trees file", sfkit::assert::light);
//     }

//     TSKitTreeSequence(tsk_treeseq_t const& tree_sequence) : _tree_sequence(tree_sequence) {}

//     ~TSKitTreeSequence() {
//         tsk_treeseq_free(&_tree_sequence);
//     }

//     TSKitTreeSequence(TSKitTreeSequence const& other)            = delete;
//     TSKitTreeSequence& operator=(TSKitTreeSequence const& other) = delete;

//     TSKitTreeSequence(TSKitTreeSequence&& other) noexcept : _tree_sequence(other._tree_sequence) {
//         other._tree_sequence = tsk_treeseq_t();
//     }

//     TSKitTreeSequence& operator=(TSKitTreeSequence&& other) {
//         if (this != &other) {
//             tsk_treeseq_free(&_tree_sequence);
//             _tree_sequence       = other._tree_sequence;
//             other._tree_sequence = tsk_treeseq_t();
//         }
//         return *this;
//     }

//     NodeId num_nodes() const {
//         return asserting_cast<NodeId>(tsk_treeseq_get_num_nodes(&_tree_sequence));
//     }

//     EdgeId num_edges() const {
//         return asserting_cast<EdgeId>(tsk_treeseq_get_num_edges(&_tree_sequence));
//     }

//     SiteId num_sites() const {
//         return asserting_cast<SiteId>(tsk_treeseq_get_num_sites(&_tree_sequence));
//     }

//     MutationId num_mutations() const {
//         return asserting_cast<MutationId>(tsk_treeseq_get_num_mutations(&_tree_sequence));
//     }

//     std::size_t num_populations() const {
//         return tsk_treeseq_get_num_populations(&_tree_sequence);
//     }

//     std::size_t num_individuals() const {
//         return tsk_treeseq_get_num_individuals(&_tree_sequence);
//     }

//     TreeId num_trees() const {
//         return asserting_cast<TreeId>(tsk_treeseq_get_num_trees(&_tree_sequence));
//     }

//     SampleId num_samples() const {
//         return asserting_cast<SampleId>(tsk_treeseq_get_num_samples(&_tree_sequence));
//     }

//     char const* file_uuid() const {
//         return tsk_treeseq_get_file_uuid(&_tree_sequence);
//     }

//     double sequence_length() const {
//         return tsk_treeseq_get_sequence_length(&_tree_sequence);
//     }

//     tsk_treeseq_t& underlying() {
//         return _tree_sequence;
//     }

//     tsk_treeseq_t const& underlying() const {
//         return _tree_sequence;
//     }

//     [[nodiscard]] bool sample_ids_are_consecutive() const {
//         // According to the tskit documentation: The array is owned by the tree sequence and should not be modified
//         // or free’d.
//         tsk_id_t const*  samples     = tsk_treeseq_get_samples(&underlying());
//         tsk_size_t const num_samples = tsk_treeseq_get_num_samples(&underlying());
//         for (tsk_id_t id = 0; id < asserting_cast<tsk_id_t>(num_samples); ++id) {
//             if (id != *(samples + id)) {
//                 return false;
//             }
//         }
//         return true;
//     }

//     bool is_sample(tsk_id_t u) const {
//         return tsk_treeseq_is_sample(&_tree_sequence, u);
//     }

//     bool is_discrete_genome() const {
//         return tsk_treeseq_get_discrete_genome(&_tree_sequence);
//     }

//     double diversity() const {
//         double                pi;
//         std::vector<tsk_id_t> samples;
//         samples.resize(num_samples());
//         std::iota(samples.begin(), samples.end(), 0);
//         tsk_size_t sample_set_sizes[] = {num_samples()};

//         auto ret =
//             tsk_treeseq_diversity(&_tree_sequence, 1, sample_set_sizes, samples.data(), 0, NULL, TSK_STAT_SITE, &pi);
//         KASSERT(ret == 0, "Failed to compute the diversity.", sfkit::assert::light);

//         return pi;
//     }

//     double diversity(SampleSet const& samples) const {
//         double     pi;
//         auto const tsk_samples        = samples.to_tsk_samples();
//         tsk_size_t sample_set_sizes[] = {tsk_samples.size()};

//         auto ret = tsk_treeseq_diversity(
//             &_tree_sequence,
//             1,
//             sample_set_sizes,
//             tsk_samples.data(),
//             0,
//             NULL,
//             TSK_STAT_SITE,
//             &pi
//         );
//         KASSERT(ret == 0, "Failed to compute the diversity.", sfkit::assert::light);

//         return pi;
//     }

//     double divergence(SampleSet const& sample_set_1, SampleSet const& sample_set_2) const {
//         double                divergence;
//         constexpr int         num_windows     = 0;
//         constexpr int         num_sample_sets = 2;
//         std::vector<tsk_id_t> samples_sets;

//         for (auto const& sample: sample_set_1) {
//             samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
//         }
//         for (auto const& sample: sample_set_2) {
//             samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
//         }
//         tsk_size_t sample_set_sizes[] = {sample_set_1.popcount(), sample_set_2.popcount()};

//         tsk_id_t      set_indexes[]    = {0, 1};
//         constexpr int num_index_tuples = 1;

//         auto ret = tsk_treeseq_divergence(
//             &_tree_sequence,
//             num_sample_sets,
//             sample_set_sizes,
//             samples_sets.data(),
//             num_index_tuples,
//             set_indexes,
//             num_windows,
//             NULL,
//             TSK_STAT_SITE,
//             &divergence
//         );
//         KASSERT(ret == 0, "Failed to compute the divergence.", sfkit::assert::light);

//         return divergence;
//     }

//     double
//     f4(SampleSet const& sample_set_1,
//        SampleSet const& sample_set_2,
//        SampleSet const& sample_set_3,
//        SampleSet const& sample_set_4) const {
//         double                f4;
//         constexpr int         num_windows     = 0;
//         constexpr int         num_sample_sets = 4;
//         std::vector<tsk_id_t> samples_sets;

//         for (auto const& sample: sample_set_1) {
//             samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
//         }
//         for (auto const& sample: sample_set_2) {
//             samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
//         }
//         for (auto const& sample: sample_set_3) {
//             samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
//         }
//         for (auto const& sample: sample_set_4) {
//             samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
//         }
//         tsk_size_t sample_set_sizes[] =
//             {sample_set_1.popcount(), sample_set_2.popcount(), sample_set_3.popcount(), sample_set_4.popcount()};

//         tsk_id_t      set_indexes[]    = {0, 1, 2, 3};
//         constexpr int num_index_tuples = 1;

//         auto ret = tsk_treeseq_f4(
//             &_tree_sequence,
//             num_sample_sets,
//             sample_set_sizes,
//             samples_sets.data(),
//             num_index_tuples,
//             set_indexes,
//             num_windows,
//             NULL,
//             TSK_STAT_SITE,
//             &f4
//         );
//         KASSERT(ret == 0, "Failed to compute the f4.", sfkit::assert::light);

//         return f4;
//     }

//     double f3(SampleSet const& sample_set_1, SampleSet const& sample_set_2, SampleSet const& sample_set_3) const {
//         double                f3;
//         constexpr int         num_windows     = 0;
//         constexpr int         num_sample_sets = 3;
//         std::vector<tsk_id_t> samples_sets;

//         for (auto const& sample: sample_set_1) {
//             samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
//         }
//         for (auto const& sample: sample_set_2) {
//             samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
//         }
//         for (auto const& sample: sample_set_3) {
//             samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
//         }
//         tsk_size_t sample_set_sizes[] = {sample_set_1.popcount(), sample_set_2.popcount(), sample_set_3.popcount()};

//         tsk_id_t      set_indexes[]    = {0, 1, 2, 3};
//         constexpr int num_index_tuples = 1;

//         auto ret = tsk_treeseq_f3(
//             &_tree_sequence,
//             num_sample_sets,
//             sample_set_sizes,
//             samples_sets.data(),
//             num_index_tuples,
//             set_indexes,
//             num_windows,
//             NULL,
//             TSK_STAT_SITE,
//             &f3
//         );
//         KASSERT(ret == 0, "Failed to compute the f3.", sfkit::assert::light);

//         return f3;
//     }

//     double f2(SampleSet const& sample_set_1, SampleSet const& sample_set_2) const {
//         double                f2;
//         constexpr int         num_windows     = 0;
//         constexpr int         num_sample_sets = 2;
//         std::vector<tsk_id_t> samples_sets;

//         for (auto const& sample: sample_set_1) {
//             samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
//         }
//         for (auto const& sample: sample_set_2) {
//             samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
//         }
//         tsk_size_t sample_set_sizes[] = {sample_set_1.popcount(), sample_set_2.popcount()};

//         tsk_id_t      set_indexes[]    = {0, 1, 2, 3};
//         constexpr int num_index_tuples = 1;

//         auto ret = tsk_treeseq_f2(
//             &_tree_sequence,
//             num_sample_sets,
//             sample_set_sizes,
//             samples_sets.data(),
//             num_index_tuples,
//             set_indexes,
//             num_windows,
//             NULL,
//             TSK_STAT_SITE,
//             &f2
//         );
//         KASSERT(ret == 0, "Failed to compute the f2.", sfkit::assert::light);

//         return f2;
//     }

//     double num_segregating_sites() const {
//         double                num_seg_sites;
//         std::vector<tsk_id_t> samples;
//         samples.resize(num_samples());

//         std::iota(samples.begin(), samples.end(), 0);
//         tsk_size_t sample_set_sizes[] = {num_samples()};

//         auto ret = tsk_treeseq_segregating_sites(
//             &_tree_sequence,
//             1,
//             sample_set_sizes,
//             samples.data(),
//             0,
//             NULL,
//             TSK_STAT_SITE,
//             &num_seg_sites
//         );
//         KASSERT(ret == 0, "Failed to compute the diversity.", sfkit::assert::light);

//         return num_seg_sites;
//     }

//     std::vector<double> allele_frequency_spectrum() const {
//         std::vector<double>   results;
//         std::vector<tsk_id_t> samples;

//         // tskit also outputs an extra 0 (for unwindowed AFS) in the first element of the results vector.
//         results.resize(num_samples() + 1);
//         samples.resize(num_samples());
//         std::iota(samples.begin(), samples.end(), 0);
//         std::array<tsk_size_t, 2> sample_set_sizes = {samples.size(), 0};

//         auto ret = tsk_treeseq_allele_frequency_spectrum(
//             &_tree_sequence,
//             1,
//             sample_set_sizes.data(),
//             samples.data(),
//             0,
//             nullptr,
//             TSK_STAT_POLARISED,
//             // reinterpret_cast<double*>(results.data())
//             results.data()
//         );
//         KASSERT(tskit_noerr(ret), "Failed to compute the allele frequency spectrum", sfkit::assert::light);
//         return results;
//     }

//     TskMutationView mutations() const {
//         return std::span(_tree_sequence.site_mutations_mem, num_mutations());
//     }

//     std::span<tsk_site_t const> sites() const {
//         KASSERT(_tree_sequence.tree_sites_mem != nullptr);
//         return std::span(_tree_sequence.tree_sites_mem, asserting_cast<size_t>(num_sites()));
//     }

//     std::span<double const> const breakpoints() const {
//         double const* breakpoints_ptr = tsk_treeseq_get_breakpoints(&_tree_sequence);
//         return std::span(breakpoints_ptr, num_trees());
//     }

//     double position_of(tsk_id_t site_id) const {
//         tsk_site_t ts_site;
//         tsk_treeseq_get_site(&_tree_sequence, site_id, &ts_site);
//         return ts_site.position;
//     }

// private:
//     bool tskit_noerr(int ret) const {
//         return ret == 0;
//     }

//     std::string   _trees_file;
//     tsk_treeseq_t _tree_sequence;
// };

// // TODO Write a proper iterator
// class TSKitTree {
// public:
//     TSKitTree(TSKitTreeSequence& tree_sequence) : _tree_sequence(tree_sequence) {
//         auto ret = tsk_tree_init(&_tree, &(_tree_sequence.underlying()), 0);
//         KASSERT(ret == 0, "Failed to initialize the tree data structure", sfkit::assert::light);

//         // Load the first tree
//         first();
//         _tree_id = 0;
//     }

//     ~TSKitTree() {
//         tsk_tree_free(&_tree);
//     }

//     bool first() {
//         _state = tsk_tree_first(&_tree);
//         KASSERT(_state >= 0, "Failed to goto the first tree.", sfkit::assert::light);
//         return _state == TSK_TREE_OK;
//     }

//     bool next() {
//         _state = tsk_tree_next(&_tree);
//         ++_tree_id;
//         KASSERT(_state >= 0, "Failed to goto the next tree.", sfkit::assert::light);
//         return _state == TSK_TREE_OK;
//     }

//     tsk_id_t tree_id() const {
//         return _tree_id;
//     }

//     bool is_valid() const {
//         return is_tree() || is_null();
//     }

//     bool is_tree() const {
//         return _state == TSK_TREE_OK;
//     }

//     bool is_null() const {
//         return _state == TSK_NULL_TREE;
//     }

//     std::size_t num_roots() const {
//         KASSERT(_state >= 0, "The tree is not valid.", sfkit::assert::light);
//         return tsk_tree_get_num_roots(&_tree);
//     }

//     std::size_t max_node_id() const {
//         return _tree_sequence.num_nodes();
//     }

//     std::size_t num_samples() const {
//         return _tree_sequence.num_samples();
//     }

//     bool is_sample(tsk_id_t node) const {
//         KASSERT(_state >= 0, "The tree is not valid.", sfkit::assert::light);
//         KASSERT(asserting_cast<tsk_size_t>(node) <= max_node_id(), "The node is not valid.", sfkit::assert::light);
//         return _tree_sequence.is_sample(node);
//     }

//     tsk_id_t num_children(tsk_id_t node) const {
//         KASSERT(_state >= 0, "The tree is not valid.", sfkit::assert::light);
//         KASSERT(asserting_cast<tsk_size_t>(node) <= max_node_id(), "The node is not valid.", sfkit::assert::light);
//         return _tree.num_children[node];
//     }

//     tsk_id_t root() const {
//         KASSERT(_state >= 0, "The tree is not valid.", sfkit::assert::light);
//         KASSERT(
//             tsk_tree_get_num_roots(&_tree) == 1ul,
//             "The tskit tree has more than one root. We do not support this yet."
//         );
//         KASSERT(_tree.virtual_root != TSK_NULL);
//         KASSERT(_tree.right_child[_tree.virtual_root] != TSK_NULL);
//         KASSERT(_tree.left_child[_tree.virtual_root] == _tree.right_child[_tree.virtual_root]);
//         KASSERT(num_children(_tree.virtual_root) == 1);
//         return _tree.right_child[_tree.virtual_root];
//     }

//     bool is_root(tsk_id_t const node) const {
//         // The following would NOT work, as if the node is not in the current tree, the parent is also set to TKS_NULL!
//         // return _tree.parent[node] == TSK_NULL;

//         return root() == node;
//     }

//     auto postorder() {
//         int        ret;
//         tsk_size_t num_nodes;
//         _postorder_nodes_resize();
//         ret = tsk_tree_postorder(&_tree, _postorder_nodes.data(), &num_nodes);
//         KASSERT(ret == 0, "Failed to get the postorder traversal of the tree.", sfkit::assert::light);
//         KASSERT(
//             num_nodes <= _postorder_nodes.size(),
//             "The number of nodes in the postorder traversal is too large to fit into the buffer.",
//             sfkit::assert::light
//         );
//         return std::span{_postorder_nodes}.subspan(0, num_nodes);
//     }

//     class EulertourView {
//     public:
//         EulertourView(tsk_tree_t const& tree, TSKitTreeSequence const& tree_sequence)
//             : _tree(tree),
//               _tree_sequence(tree_sequence) {}
//         class iterator {
//         public:
//             using iterator_category = std::forward_iterator_tag;
//             using difference_type   = std::ptrdiff_t;
//             using value_type        = tsk_id_t;
//             using pointer           = value_type*;
//             using reference         = value_type&;
//             struct sentinel {};

//             iterator(tsk_tree_t const& tree, TSKitTreeSequence const& tree_sequence)
//                 : _tree(tree),
//                   _just_moved_up(false),
//                   _tree_sequence(tree_sequence) {
//                 KASSERT(
//                     _tree.virtual_root != TSK_NULL,
//                     "The virtual root of the tree does not exist.",
//                     sfkit::assert::light
//                 );
//                 KASSERT(
//                     _tree.right_child[_tree.virtual_root] != TSK_NULL,
//                     "The virtual root of the tree has no children.",
//                     sfkit::assert::light
//                 );
//                 KASSERT(
//                     _tree.left_child[_tree.virtual_root] == _tree.right_child[_tree.virtual_root],
//                     "The virtual root of the tree has more than one child.",
//                     sfkit::assert::light
//                 );
//                 _current_node = _tree.left_child[_tree.virtual_root];
//             }

//             iterator& operator++() {
//                 if (_just_moved_up || is_sample()) {
//                     // The children of this node were already processed OR there are no children because the current
//                     // node is a sample node.
//                     if (_tree.right_sib[_current_node] != TSK_NULL) {
//                         // Move to the next unprocessed sibling of this node.
//                         _current_node  = _tree.right_sib[_current_node];
//                         _just_moved_up = false;
//                     } else {
//                         // There are no more siblings of this node -> move up.
//                         _current_node  = _tree.parent[_current_node];
//                         _just_moved_up = true;
//                     }
//                 } else {
//                     // We just moved down or right (to a sibling), the children of this node are not processed yet
//                     //   -> move down
//                     _current_node = _tree.left_child[_current_node];
//                 }
//                 return *this;
//             }

//             iterator operator++(int) {
//                 iterator tmp(*this);

//                 this->operator++();
//                 return tmp;
//             }

//             [[nodiscard]] bool operator==(iterator const& other) const {
//                 return _current_node == other._current_node && &_tree == &other._tree
//                        && _just_moved_up == other._just_moved_up;
//             }

//             [[nodiscard]] bool operator!=(iterator const& other) const {
//                 return !(*this == other);
//             }

//             [[nodiscard]] bool operator==(sentinel) {
//                 return _current_node == TSK_NULL;
//             }

//             [[nodiscard]] reference operator*() {
//                 KASSERT(_current_node != TSK_NULL);
//                 return _current_node;
//             }

//             [[nodiscard]] reference node_id() {
//                 return operator*();
//             }

//             [[nodiscard]] bool first_visit() {
//                 return !_just_moved_up;
//             }

//             [[nodiscard]] bool second_visit() {
//                 return _just_moved_up;
//             }

//             [[nodiscard]] bool is_sample() {
//                 KASSERT(_current_node != TSK_NULL);
//                 return _tree_sequence.is_sample(_current_node);
//             }

//             [[nodiscard]] pointer operator->() {
//                 KASSERT(_current_node != TSK_NULL);
//                 return &_current_node;
//             }

//         private:
//             tsk_tree_t const&        _tree;
//             tsk_id_t                 _current_node;
//             bool                     _just_moved_up;
//             TSKitTreeSequence const& _tree_sequence;
//         };

//         auto begin() const {
//             return iterator{_tree, _tree_sequence};
//         }

//         auto end() const {
//             return iterator::sentinel{};
//         }

//     private:
//         tsk_tree_t const&        _tree;
//         TSKitTreeSequence const& _tree_sequence;
//     };

//     EulertourView eulertour() const {
//         return EulertourView{_tree, _tree_sequence};
//     }

//     class Children {
//     public:
//         Children(tsk_tree_t const& tree, tsk_id_t const parent) : _tree(tree), _parent(parent) {}
//         class children_iterator {
//         public:
//             using iterator_category = std::forward_iterator_tag;
//             using difference_type   = std::ptrdiff_t;
//             using value_type        = tsk_id_t;
//             using pointer           = value_type*;
//             using reference         = value_type&;
//             struct sentinel {};

//             children_iterator(tsk_tree_t const& tree, tsk_id_t const parent)
//                 : _tree(tree),
//                   _child(tree.left_child[parent]) {}

//             children_iterator& operator++() {
//                 _child = _tree.right_sib[_child];
//                 // This provided a small but measurable speedup (~10%)
//                 __builtin_prefetch(&_tree.right_sib[_child]);
//                 return *this;
//             }

//             children_iterator operator++(int) {
//                 children_iterator tmp(*this);
//                 this->            operator++();
//                 return tmp;
//             }

//             bool operator==(children_iterator const& other) const {
//                 return _child == other._child && &_tree == &other._tree;
//             }

//             bool operator!=(children_iterator const& other) const {
//                 return !(*this == other);
//             }

//             bool operator==(sentinel) {
//                 return _child == TSK_NULL;
//             }

//             reference operator*() {
//                 return _child;
//             }

//             pointer operator->() {
//                 return &_child;
//             }

//         private:
//             tsk_tree_t const& _tree;
//             tsk_id_t          _child;
//         };

//         auto begin() const {
//             return children_iterator{_tree, _parent};
//         }

//         auto end() const {
//             return children_iterator::sentinel{};
//         }

//     private:
//         tsk_tree_t const& _tree;
//         tsk_id_t const    _parent;
//     };

//     Children children(tsk_id_t const parent) {
//         return Children{_tree, parent};
//     }

//     tsk_id_t parent(tsk_id_t const node) const {
//         KASSERT(_state >= 0, "The tree is not valid.", sfkit::assert::light);
//         KASSERT(asserting_cast<tsk_size_t>(node) <= max_node_id(), "The node is not valid.", sfkit::assert::light);
//         return _tree.parent[node];
//     }

// private:
//     size_t _current_tree_size_bound() const {
//         return tsk_tree_get_size_bound(&_tree);
//     }

//     // Reserve space for the traversal lists. Sadly, tsk provides us only an upper bound /on this tree/.
//     // We therefore might need to resize the vectors later on.
//     void _postorder_nodes_resize() {
//         if (_postorder_nodes.size() < _current_tree_size_bound()) { // We never want to shrink the vector
//             _postorder_nodes.resize(_current_tree_size_bound(), TSK_NULL);
//         }
//     }

//     // TODO Create a variadic templated _crate_tsk_sample_sets and use it in the functions above.

//     TSKitTreeSequence&    _tree_sequence;
//     tsk_tree_t            _tree;
//     tsk_id_t              _tree_id;
//     int                   _state = -1;
//     std::vector<tsk_id_t> _postorder_nodes;
// };
// } // namespace sfkit::tskit

#include "ChangedNodesView.hpp"
#include "ChildrenView.hpp"
#include "EulertourView.hpp"
#include "TSKitTree.hpp"
#include "TSKitTreeSequence.hpp"
