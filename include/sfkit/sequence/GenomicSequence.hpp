#pragma once

#include <algorithm>
#include <cstddef>
#include <span>
#include <unordered_set>
#include <vector>

#include <kassert/kassert.hpp>
#include <tskit/core.h>
#include <tskit/tables.h>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/graph/types.hpp"
#include "sfkit/io/vector_serialization.hpp"
#include "sfkit/sequence/Mutation.hpp"
#include "sfkit/sequence/Sequence.hpp"
#include "sfkit/sequence/TSKitSiteToTreeMapper.hpp"
#include "sfkit/tskit/tskit.hpp"

namespace sfkit::sequence {

class GenomicSequence {
public:
    GenomicSequence(SiteId num_sites_hint = 0, MutationId num_mutations_hint = 0) {
        _sites.reserve(asserting_cast<size_t>(num_sites_hint));
        _mutation_indices.reserve(asserting_cast<size_t>(num_sites_hint));
        _mutations.reserve(num_mutations_hint);
    }

    [[nodiscard]] SiteId num_sites() const {
        return asserting_cast<SiteId>(_sites.size());
    }

    [[nodiscard]] MutationId num_mutations() const {
        return asserting_cast<MutationId>(_mutations.size());
    }

    [[nodiscard]] AllelicState ancestral_state(SiteId site_id) const {
        KASSERT(site_id >= 0, "Site ID is invalid.", sfkit::assert::light);
        KASSERT(asserting_cast<size_t>(site_id) < _sites.size(), "Site ID is out of bounds", sfkit::assert::light);
        return _sites[asserting_cast<size_t>(site_id)];
    }

    [[nodiscard]] AllelicState& ancestral_state(SiteId site_id) {
        KASSERT(site_id >= 0, "Site ID is invalid.", sfkit::assert::light);
        KASSERT(asserting_cast<size_t>(site_id) < _sites.size(), "Site ID is out of bounds", sfkit::assert::light);
        return _sites[asserting_cast<size_t>(site_id)];
    }

    [[nodiscard]] Mutation const& mutation_by_id(size_t mutation_id) const {
        KASSERT(mutation_id < _mutations.size(), "Mutation ID is out of bounds", sfkit::assert::light);
        return _mutations[mutation_id];
    }

    void set(SiteId site_id, AllelicState state) {
        KASSERT(site_id >= 0, "Site ID is invalid.", sfkit::assert::light);
        KASSERT(asserting_cast<size_t>(site_id) < _sites.size(), "Site ID is out of bounds", sfkit::assert::light);
        _sites[asserting_cast<size_t>(site_id)] = state;
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
        _mutation_indices_valid = false;
    }

    void emplace_back(Mutation&& mutation) {
        _mutations.emplace_back(std::move(mutation));
        _mutation_indices_valid = false;
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
        KASSERT(site_id >= 0, "Site ID is invalid.", sfkit::assert::light);
        return _sites[asserting_cast<size_t>(site_id)];
    }

    void build_mutation_indices() {
        KASSERT(mutations_are_sorted_by_site(), "Mutations are not sorted by site.", sfkit::assert::light);
        _mutation_indices.clear();
        _mutation_indices.reserve(_sites.size());
        MutationId mutation_idx = 0;
        _mutation_indices.push_back(0);
        for (SiteId site_id = 0; site_id < num_sites(); ++site_id) {
            while (mutation_idx < _mutations.size()
                   && _mutations[asserting_cast<size_t>(mutation_idx)].site_id() == site_id) {
                ++mutation_idx;
            }
            _mutation_indices.push_back(mutation_idx);
        }
        KASSERT(
            _mutations.size() == asserting_cast<size_t>(mutation_idx),
            "Mutation index is not at the end of the mutation vector",
            sfkit::assert::light
        );
        // Add sentinel
        _mutation_indices.push_back(mutation_idx);
        _mutation_indices_valid = true;
    }

    [[nodiscard]] bool mutation_indices_are_built() const {
        return _mutation_indices_valid;
    }

    [[nodiscard]] bool mutations_are_sorted_by_site() const {
        return std::is_sorted(_mutations.begin(), _mutations.end(), [](Mutation const& lhs, Mutation const& rhs) {
            return lhs.site_id() < rhs.site_id();
        });
    }

    [[nodiscard]] bool mutations_are_sorted_by_tree_id() const {
        return std::is_sorted(_mutations.begin(), _mutations.end(), [](Mutation const& lhs, Mutation const& rhs) {
            // TODO Rename subtree_id() to tree_id()
            return lhs.node_id() < rhs.node_id();
        });
    }

    [[nodiscard]] MutationView mutations_at_site(SiteId const site_id) const {
        KASSERT(site_id >= 0, "Site ID is invalid.", sfkit::assert::light);
        KASSERT(_mutation_indices_valid, "Mutations indices need to be rebuild first.", sfkit::assert::light);
        KASSERT(
            asserting_cast<size_t>(site_id) < _mutation_indices.size() + 1,
            "Site ID is out of bounds",
            sfkit::assert::light
        );
        KASSERT(
            _mutations.size() == _mutation_indices.back(),
            "The _mutations_indices sentinel seems to be broken.",
            sfkit::assert::light
        );
        return std::span(_mutations)
            .subspan(
                _mutation_indices[asserting_cast<size_t>(site_id)],
                _mutation_indices[asserting_cast<size_t>(site_id + 1)]
                    - _mutation_indices[asserting_cast<size_t>(site_id)]
            );
    }

    [[nodiscard]] std::unordered_set<NodeId> subtrees_with_mutations() const {
        std::unordered_set<NodeId> subtrees_with_mutations;

        std::for_each(_mutations.begin(), _mutations.end(), [&subtrees_with_mutations](Mutation const& mutation) {
            subtrees_with_mutations.insert(mutation.node_id());
        });

        return subtrees_with_mutations;
    }

    template <class Archive>
    void serialize(Archive& archive) {
        build_mutation_indices();
        archive(_sites, _mutation_indices, _mutation_indices_valid, _mutations);
    }

    void save(std::ostream& os) {
        build_mutation_indices();
        sfkit::io::utils::serialize(os, _sites);
        sfkit::io::utils::serialize(os, _mutation_indices);
        sfkit::io::utils::serialize(os, _mutations);
        os.write(reinterpret_cast<char const*>(&_mutation_indices_valid), sizeof(_mutation_indices_valid));
    }

    void load(std::istream& is) {
        sfkit::io::utils::deserialize(is, _sites);
        sfkit::io::utils::deserialize(is, _mutation_indices);
        sfkit::io::utils::deserialize(is, _mutations);
        is.read(reinterpret_cast<char*>(&_mutation_indices_valid), sizeof(_mutation_indices_valid));
    }

private:
    std::vector<AllelicState> _sites;
    std::vector<MutationId>   _mutation_indices; // Maps SiteId to MutationId
    std::vector<Mutation>     _mutations;
    bool                      _mutation_indices_valid = false;
};
} // namespace sfkit::sequence

inline std::ostream& operator<<(std::ostream& os, sfkit::sequence::MutationView const& mutations) {
    os << "{ ";
    for (auto const& mutation: mutations) {
        os << "Mutation<site:" << mutation.site_id() << ", derived:" << mutation.allelic_state() << "> ";
    }
    os << " }";
    return os;
}

inline std::ostream& operator<<(std::ostream& os, sfkit::sequence::Mutation const& mutation) {
    os << "Mutation<site:" << mutation.site_id() << ", derived:" << mutation.allelic_state() << ">";
    return os;
}
