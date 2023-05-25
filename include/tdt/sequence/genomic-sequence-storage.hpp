#pragma once

#include <algorithm>
#include <cstddef>
#include <span>
#include <vector>

#include <kassert/kassert.hpp>
#include <tskit/tables.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/graph/common.hpp"
#include "tdt/load/forest-compressor.hpp"
#include "tdt/sequence/tskit-site-to-tree-mapper.hpp"
#include "tdt/tskit.hpp"

using SiteId = size_t;
// TODO Use only 2 bits per site
using AllelicState                         = char;
constexpr AllelicState InvalidAllelicState = -1;

// TODO Implement something safer here (e.g. what happens if '1' and 'A' gets mixed in the same dataset?)
class PerfectNumericHasher {
public:
    static inline unsigned char to_idx(AllelicState state) {
        char const idx = state - '0';
        KASSERT((idx >= 0 && idx < num_states), "Invalid state: " << state, tdt::assert::light);
        return static_cast<unsigned char>(idx);
    }

    static constexpr unsigned int num_states = 4;
};

class PerfectDNAHasher {
public:
    static inline unsigned char to_idx(AllelicState state) {
        // The characters A, C, T and G differ pairwise in the second and third least significant bit.
        // A 0100 0001
        // C 0100 0011
        // G 0100 0111
        // T 0101 0100
        //         ^^
        constexpr char          ISOLATE_BITMASK = 0b0000'0110;
        constexpr unsigned char SHIFT           = 1;
        unsigned char           idx             = asserting_cast<unsigned char>(state & ISOLATE_BITMASK) >> SHIFT;
        KASSERT(idx < num_states, "Invalid state: " << state, tdt::assert::light);
        return idx;
    }

    static constexpr unsigned int num_states = 4;
};

class GenomicSequenceStorage {
public:
    class Mutation {
    public:
        Mutation(
            SiteId       site_id,
            TreeId       tree_id,
            SubtreeId    subtree_id,
            AllelicState state,
            tsk_id_t     parent_mutation_id = TSK_NULL
        ) noexcept
            : _site_id(site_id),
              _derived_state(state),
              _tree_id(tree_id),
              _subtree_id(subtree_id),
              _parent_mutation_id(parent_mutation_id) {}

        [[nodiscard]] SiteId site_id() const {
            return _site_id;
        }

        [[nodiscard]] AllelicState allelic_state() const {
            return _derived_state;
        }

        [[nodiscard]] SubtreeId tree_id() const {
            return _tree_id;
        }

        [[nodiscard]] TreeId subtree_id() const {
            return _subtree_id;
        }

        // TODO Use a reference instead of an id?
        [[nodiscard]] tsk_id_t parent_mutation_id() const {
            return _parent_mutation_id;
        }

    private:
        // TODO Compress this representation. Through the use of indices we could get rid of the tree id
        // maybe even the site id, using a retrieval data structure? We also know the mapping site_id -> tree_id
        SiteId       _site_id;
        AllelicState _derived_state;
        TreeId       _tree_id;
        SubtreeId    _subtree_id; // TODO Change to compressed forest node id
        tsk_id_t     _parent_mutation_id;
    };

    GenomicSequenceStorage(size_t num_sites_hint = 0, size_t num_mutations_hint = 0) {
        _sites.reserve(num_sites_hint);
        _mutation_indices.reserve(num_sites_hint);
        _mutations.reserve(num_mutations_hint);
    }

    GenomicSequenceStorage(TSKitTreeSequence const& tree_sequence, CompressedForest const& compressed_forest)
        : _num_sites(tree_sequence.num_sites()) {
        TSKitSiteToTreeMapper site2tree(tree_sequence);

        // Store ancestral states
        for (auto&& site: tree_sequence.sites()) {
            KASSERT(site.ancestral_state_length == 1u, "Ancestral state length is not 1", tdt::assert::light);
            _sites.emplace_back(*site.ancestral_state);
        }
        KASSERT(
            _sites.size() == tree_sequence.num_sites(),
            "Number of sites reported by num_sites() and in the sites() iterator does not match",
            tdt::assert::light
        );

        // Store mutations
        KASSERT(
            std::is_sorted(
                tree_sequence.mutations().begin(),
                tree_sequence.mutations().end(),
                [](auto&& lhs, auto&& rhs) { return lhs.site < rhs.site; }
            ),
            "Mutations are not sorted by site",
            tdt::assert::normal
        );

        // The mutations are sorted by site.
        for (auto& ts_mutation: tree_sequence.mutations()) {
            // TODO Use less bits for site and tree ids
            SiteId const    site_id       = asserting_cast<SiteId>(ts_mutation.site);
            TreeId const    tree_id       = site2tree(ts_mutation.site);
            SubtreeId const cf_subtree_id = compressed_forest.ts_node2cf_subtree(tree_id, ts_mutation.node);

            KASSERT(ts_mutation.derived_state_length == 1u, "Derived state length is not 1", tdt::assert::light);
            AllelicState const derived_state = *ts_mutation.derived_state;

            KASSERT(ts_mutation.node != TSK_NULL, "Mutation node is null", tdt::assert::light);

            tsk_id_t parent_mutation_id = ts_mutation.parent;

            tsk_id_t mutation_id = asserting_cast<tsk_id_t>(_mutations.size());
            KASSERT(
                mutation_id == ts_mutation.id,
                "Mutation ID is not equal to the index in the mutations vector",
                tdt::assert::light
            );

            _mutations.emplace_back(site_id, tree_id, cf_subtree_id, derived_state, parent_mutation_id);
        }

        build_mutation_indices();
    }

    [[nodiscard]] size_t num_sites() const {
        return _sites.size();
    }

    [[nodiscard]] size_t num_mutations() const {
        return _mutations.size();
    }

    [[nodiscard]] AllelicState ancestral_state(SiteId site_id) const {
        KASSERT(site_id < _sites.size(), "Site ID is out of bounds", tdt::assert::light);
        return _sites[site_id];
    }

    [[nodiscard]] AllelicState& ancestral_state(SiteId site_id) {
        KASSERT(site_id < _sites.size(), "Site ID is out of bounds", tdt::assert::light);
        return _sites[site_id];
    }

    [[nodiscard]] Mutation const& mutation_by_id(size_t mutation_id) const {
        KASSERT(mutation_id < _mutations.size(), "Mutation ID is out of bounds", tdt::assert::light);
        return _mutations[mutation_id];
    }

    void set(SiteId site_id, AllelicState state) {
        KASSERT(site_id < _sites.size(), "Site ID is out of bounds", tdt::assert::light);
        _sites[site_id] = state;
    }

    void push_back(AllelicState state) {
        _sites.push_back(state);
    }

    void emplace_back(AllelicState&& state) {
        _sites.emplace_back(std::move(state));
    }

    template <class... Args>
    requires std::constructible_from<SiteId, Args...>
    void emplace_back(Args&&... args) {
        _sites.emplace_back(std::forward<Args>(args)...);
    }

    void push_back(Mutation const& mutation) {
        _mutations.push_back(mutation);
    }

    void emplace_back(Mutation&& mutation) {
        _mutations.emplace_back(std::move(mutation));
    }

    template <class... Args>
    requires std::constructible_from<Mutation, Args...>
    void emplace_back(Args&&... args) {
        _mutations.emplace_back(std::forward<Args>(args)...);
    }

    [[nodiscard]] AllelicState operator[](SiteId site_id) const {
        return ancestral_state(site_id);
    }

    [[nodiscard]] AllelicState& operator[](SiteId site_id) {
        return _sites[site_id];
    }

    void build_mutation_indices() {
        KASSERT(mutations_are_sorted_by_site(), "Mutations are not sorted by site.", tdt::assert::light);
        _mutation_indices.clear();
        _mutation_indices.reserve(_sites.size());
        // TODO What if the sites are not set  / are set?
        size_t mutation_idx = 0;
        _mutation_indices.push_back(0);
        for (SiteId site_id = 0; site_id < _num_sites; ++site_id) {
            while (mutation_idx < _mutations.size() && _mutations[mutation_idx].site_id() == site_id) {
                ++mutation_idx;
            }
            _mutation_indices.push_back(mutation_idx);
        }
        KASSERT(
            _mutations.size() == mutation_idx,
            "Mutation index is not at the end of the mutation vector",
            tdt::assert::light
        );
        // Add sentinel
        _mutation_indices.push_back(mutation_idx);
    }

    [[nodiscard]] bool mutations_are_sorted_by_site() const {
        return std::is_sorted(_mutations.begin(), _mutations.end(), [](Mutation const& lhs, Mutation const& rhs) {
            return lhs.site_id() < rhs.site_id();
        });
    }

    [[nodiscard]] bool mutations_are_sorted_by_tree_id() const {
        return std::is_sorted(_mutations.begin(), _mutations.end(), [](Mutation const& lhs, Mutation const& rhs) {
            // TODO Rename subtree_id() to tree_id()
            return lhs.subtree_id() < rhs.subtree_id();
        });
    }

    [[nodiscard]] std::span<const Mutation> mutations_at_site(SiteId const site_id) const {
        return std::span(_mutations)
            .subspan(_mutation_indices[site_id], _mutation_indices[site_id + 1] - _mutation_indices[site_id]);
    }

private:
    std::vector<AllelicState> _sites;
    std::vector<size_t>       _mutation_indices; // Maps SiteId to MutationId
    std::vector<Mutation>     _mutations;
    size_t                    _num_sites;
};

inline std::ostream& operator<<(std::ostream& os, std::span<GenomicSequenceStorage::Mutation> const& mutations) {
    os << "{ ";
    for (auto const& mutation: mutations) {
        os << "Mutation<site:" << mutation.site_id() << ", derived:" << mutation.allelic_state() << "> ";
    }
    os << " }";
    return os;
}

inline std::ostream& operator<<(std::ostream& os, GenomicSequenceStorage::Mutation const& mutation) {
    os << "Mutation<site:" << mutation.site_id() << ", derived:" << mutation.allelic_state() << ">";
    return os;
}
