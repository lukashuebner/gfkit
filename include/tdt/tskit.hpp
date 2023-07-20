#pragma once

#include <cstddef>
#include <iterator>
#include <numeric>
#include <span>

#include <kassert/kassert.hpp>
#include <tskit.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/checking_casts.hpp"
#include "tdt/graph/common.hpp"
#include "tdt/samples/sample-set.hpp"

using TskMutationView = std::span<tsk_mutation_t const>;

class TSKitTreeSequence {
public:
    TSKitTreeSequence(std::string const& trees_file) : _trees_file(trees_file) {
        int ret;

        // Load the tree sequence from the .trees file
        ret = tsk_treeseq_load(&_tree_sequence, _trees_file.c_str(), 0);
        KASSERT(tskit_noerr(ret), "Failed to load tree sequence from the .trees file", tdt::assert::light);
    }

    TSKitTreeSequence(tsk_treeseq_t const& tree_sequence) : _tree_sequence(tree_sequence) {}

    ~TSKitTreeSequence() {
        tsk_treeseq_free(&_tree_sequence);
    }

    TSKitTreeSequence(TSKitTreeSequence const& other)            = delete;
    TSKitTreeSequence& operator=(TSKitTreeSequence const& other) = delete;

    TSKitTreeSequence(TSKitTreeSequence&& other) noexcept : _tree_sequence(other._tree_sequence) {
        other._tree_sequence = tsk_treeseq_t();
    }

    TSKitTreeSequence& operator=(TSKitTreeSequence&& other) {
        if (this != &other) {
            tsk_treeseq_free(&_tree_sequence);
            _tree_sequence       = other._tree_sequence;
            other._tree_sequence = tsk_treeseq_t();
        }
        return *this;
    }

    std::size_t num_nodes() const {
        return tsk_treeseq_get_num_nodes(&_tree_sequence);
    }

    std::size_t num_edges() const {
        return tsk_treeseq_get_num_edges(&_tree_sequence);
    }

    std::size_t num_sites() const {
        return tsk_treeseq_get_num_sites(&_tree_sequence);
    }

    std::size_t num_mutations() const {
        return tsk_treeseq_get_num_mutations(&_tree_sequence);
    }

    std::size_t num_populations() const {
        return tsk_treeseq_get_num_populations(&_tree_sequence);
    }

    std::size_t num_individuals() const {
        return tsk_treeseq_get_num_individuals(&_tree_sequence);
    }

    std::size_t num_trees() const {
        return tsk_treeseq_get_num_trees(&_tree_sequence);
    }

    std::size_t num_samples() const {
        return tsk_treeseq_get_num_samples(&_tree_sequence);
    }

    char const* file_uuid() const {
        return tsk_treeseq_get_file_uuid(&_tree_sequence);
    }

    double sequence_length() const {
        return tsk_treeseq_get_sequence_length(&_tree_sequence);
    }

    tsk_treeseq_t& underlying() {
        return _tree_sequence;
    }

    tsk_treeseq_t const& underlying() const {
        return _tree_sequence;
    }

    [[nodiscard]] bool sample_ids_are_consecutive() const {
        // According to the tskit documentation: The array is owned by the tree sequence and should not be modified
        // or free’d.
        tsk_id_t const*  samples     = tsk_treeseq_get_samples(&underlying());
        tsk_size_t const num_samples = tsk_treeseq_get_num_samples(&underlying());
        for (tsk_id_t id = 0; id < asserting_cast<tsk_id_t>(num_samples); ++id) {
            if (id != *(samples + id)) {
                return false;
            }
        }
        return true;
    }

    bool is_sample(tsk_id_t u) const {
        return tsk_treeseq_is_sample(&_tree_sequence, u);
    }

    bool is_discrete_genome() const {
        return tsk_treeseq_get_discrete_genome(&_tree_sequence);
    }

    double diversity() const {
        double                pi;
        std::vector<tsk_id_t> samples;
        samples.resize(num_samples());
        std::iota(samples.begin(), samples.end(), 0);
        tsk_size_t sample_set_sizes[] = {num_samples()};

        auto ret =
            tsk_treeseq_diversity(&_tree_sequence, 1, sample_set_sizes, samples.data(), 0, NULL, TSK_STAT_SITE, &pi);
        KASSERT(ret == 0, "Failed to compute the diversity.", tdt::assert::light);

        return pi;
    }

    double diversity(SampleSet samples) const {
        double     pi;
        auto const tsk_samples        = samples.to_tsk_samples();
        tsk_size_t sample_set_sizes[] = {tsk_samples.size()};

        auto ret = tsk_treeseq_diversity(
            &_tree_sequence,
            1,
            sample_set_sizes,
            tsk_samples.data(),
            0,
            NULL,
            TSK_STAT_SITE,
            &pi
        );
        KASSERT(ret == 0, "Failed to compute the diversity.", tdt::assert::light);

        return pi;
    }

    double divergence(SampleSet const& sample_set_1, SampleSet const& sample_set_2) const {
        double                divergence;
        constexpr int         num_windows     = 0;
        constexpr int         num_sample_sets = 2;
        std::vector<tsk_id_t> samples_sets;

        for (auto const& sample: sample_set_1) {
            samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
        }
        for (auto const& sample: sample_set_2) {
            samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
        }
        tsk_size_t sample_set_sizes[] = {sample_set_1.popcount(), sample_set_2.popcount()};

        tsk_id_t      set_indexes[]    = {0, 1};
        constexpr int num_index_tuples = 1;

        auto ret = tsk_treeseq_divergence(
            &_tree_sequence,
            num_sample_sets,
            sample_set_sizes,
            samples_sets.data(),
            num_index_tuples,
            set_indexes,
            num_windows,
            NULL,
            TSK_STAT_SITE,
            &divergence
        );
        KASSERT(ret == 0, "Failed to compute the divergence.", tdt::assert::light);

        return divergence;
    }

    double num_segregating_sites() const {
        double                num_seg_sites;
        std::vector<tsk_id_t> samples;
        samples.resize(num_samples());

        std::iota(samples.begin(), samples.end(), 0);
        tsk_size_t sample_set_sizes[] = {num_samples()};

        auto ret = tsk_treeseq_segregating_sites(
            &_tree_sequence,
            1,
            sample_set_sizes,
            samples.data(),
            0,
            NULL,
            TSK_STAT_SITE,
            &num_seg_sites
        );
        KASSERT(ret == 0, "Failed to compute the diversity.", tdt::assert::light);

        return num_seg_sites;
    }

    std::vector<double> allele_frequency_spectrum() {
        std::vector<double>   results;
        std::vector<tsk_id_t> samples;

        // tskit also outputs an extra 0 (for unwindowed AFS) in the first element of the results vector.
        results.resize(num_samples() + 1);
        samples.resize(num_samples());
        std::iota(samples.begin(), samples.end(), 0);
        std::array<tsk_size_t, 2> sample_set_sizes = {samples.size(), 0};

        auto ret = tsk_treeseq_allele_frequency_spectrum(
            &_tree_sequence,
            1,
            sample_set_sizes.data(),
            samples.data(),
            0,
            nullptr,
            TSK_STAT_POLARISED,
            // reinterpret_cast<double*>(results.data())
            results.data()
        );
        KASSERT(tskit_noerr(ret), "Failed to compute the allele frequency spectrum", tdt::assert::light);
        return results;
    }

    TskMutationView mutations() const {
        return std::span(_tree_sequence.site_mutations_mem, num_mutations());
    }

    std::span<tsk_site_t const> sites() const {
        KASSERT(_tree_sequence.tree_sites_mem != nullptr);
        return std::span(_tree_sequence.tree_sites_mem, num_sites());
    }

    std::span<double const> const breakpoints() const {
        double const* breakpoints_ptr = tsk_treeseq_get_breakpoints(&_tree_sequence);
        return std::span(breakpoints_ptr, num_trees());
    }

    double position_of(tsk_id_t site_id) const {
        tsk_site_t ts_site;
        tsk_treeseq_get_site(&_tree_sequence, site_id, &ts_site);
        return ts_site.position;
    }

private:
    bool tskit_noerr(int ret) const {
        return ret == 0;
    }

    std::string   _trees_file;
    tsk_treeseq_t _tree_sequence;
};

// TODO Write a proper iterator
class TSKitTree {
public:
    TSKitTree(TSKitTreeSequence& tree_sequence) : _tree_sequence(tree_sequence) {
        auto ret = tsk_tree_init(&_tree, &(_tree_sequence.underlying()), 0);
        KASSERT(ret == 0, "Failed to initialize the tree data structure", tdt::assert::light);

        // Load the first tree
        first();
        _tree_id = 0;
    }

    ~TSKitTree() {
        tsk_tree_free(&_tree);
    }

    bool first() {
        _state = tsk_tree_first(&_tree);
        KASSERT(_state >= 0, "Failed to goto the first tree.", tdt::assert::light);
        return _state == TSK_TREE_OK;
    }

    bool next() {
        _state = tsk_tree_next(&_tree);
        ++_tree_id;
        KASSERT(_state >= 0, "Failed to goto the next tree.", tdt::assert::light);
        return _state == TSK_TREE_OK;
    }

    tsk_id_t tree_id() const {
        return _tree_id;
    }

    bool is_valid() const {
        KASSERT(_state >= 0, "The tree is not valid.", tdt::assert::light);
        return _state == TSK_TREE_OK;
    }

    std::size_t num_roots() const {
        KASSERT(_state >= 0, "The tree is not valid.", tdt::assert::light);
        return tsk_tree_get_num_roots(&_tree);
    }

    std::size_t max_node_id() const {
        return _tree_sequence.num_nodes();
    }

    std::size_t num_samples() const {
        return _tree_sequence.num_samples();
    }

    bool is_sample(tsk_id_t node) const {
        KASSERT(_state >= 0, "The tree is not valid.", tdt::assert::light);
        KASSERT(asserting_cast<tsk_size_t>(node) <= max_node_id(), "The node is not valid.", tdt::assert::light);
        return _tree_sequence.is_sample(node);
    }

    tsk_id_t num_children(tsk_id_t node) const {
        KASSERT(_state >= 0, "The tree is not valid.", tdt::assert::light);
        KASSERT(asserting_cast<tsk_size_t>(node) <= max_node_id(), "The node is not valid.", tdt::assert::light);
        return _tree.num_children[node];
    }

    tsk_id_t root() const {
        KASSERT(_state >= 0, "The tree is not valid.", tdt::assert::light);
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

    bool is_root(tsk_id_t const node) const {
        // The following would NOT work, as if the node is not in the current tree, the parent is also set to TKS_NULL!
        // return _tree.parent[node] == TSK_NULL;

        return root() == node;
    }

    auto postorder() {
        int        ret;
        tsk_size_t num_nodes;
        _postorder_nodes_resize();
        ret = tsk_tree_postorder(&_tree, _postorder_nodes.data(), &num_nodes);
        KASSERT(ret == 0, "Failed to get the postorder traversal of the tree.", tdt::assert::light);
        KASSERT(
            num_nodes <= _postorder_nodes.size(),
            "The number of nodes in the postorder traversal is too large to fit into the buffer.",
            tdt::assert::light
        );
        return std::span{_postorder_nodes}.subspan(0, num_nodes);
    }

    class Children {
    public:
        Children(tsk_tree_t const& tree, tsk_id_t const parent) : _tree(tree), _parent(parent) {}
        class children_iterator {
        public:
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using value_type        = tsk_id_t;
            using pointer           = value_type*;
            using reference         = value_type&;
            struct sentinel {};

            children_iterator(tsk_tree_t const& tree, tsk_id_t const parent)
                : _tree(tree),
                  _child(tree.left_child[parent]) {}

            children_iterator& operator++() {
                _child = _tree.right_sib[_child];
                // This provided a small but measurable speedup (~10%)
                __builtin_prefetch(&_tree.right_sib[_child]);
                return *this;
            }

            children_iterator operator++(int) {
                children_iterator tmp(*this);
                this->            operator++();
                return tmp;
            }

            bool operator==(children_iterator const& other) const {
                return _child == other._child && &_tree == &other._tree;
            }

            bool operator!=(children_iterator const& other) const {
                return !(*this == other);
            }

            bool operator==(sentinel) {
                return _child == TSK_NULL;
            }

            reference operator*() {
                return _child;
            }

            pointer operator->() {
                return &_child;
            }

        private:
            tsk_tree_t const& _tree;
            tsk_id_t          _child;
        };

        auto begin() const {
            return children_iterator{_tree, _parent};
        }

        auto end() const {
            return children_iterator::sentinel{};
        }

    private:
        tsk_tree_t const& _tree;
        tsk_id_t const    _parent;
    };

    Children children(tsk_id_t const parent) {
        return Children{_tree, parent};
    }

    tsk_id_t parent(tsk_id_t const node) const {
        KASSERT(_state >= 0, "The tree is not valid.", tdt::assert::light);
        KASSERT(asserting_cast<tsk_size_t>(node) <= max_node_id(), "The node is not valid.", tdt::assert::light);
        return _tree.parent[node];
    }

private:
    size_t _current_tree_size_bound() const {
        return tsk_tree_get_size_bound(&_tree);
    }

    // Reserve space for the traversal lists. Sadly, tsk provides us only an upper bound /on this tree/.
    // We therefore might need to resize the vectors later on.
    void _postorder_nodes_resize() {
        if (_postorder_nodes.size() < _current_tree_size_bound()) { // We never want to shrink the vector
            _postorder_nodes.resize(_current_tree_size_bound(), TSK_NULL);
        }
    }

    TSKitTreeSequence&    _tree_sequence;
    tsk_tree_t            _tree;
    tsk_id_t              _tree_id;
    int                   _state = -1;
    std::vector<tsk_id_t> _postorder_nodes;
};
