#pragma once

#include "fmt/core.h"
#include <cstddef>
#include <vector>

#include <kassert/kassert.hpp>

#include "tdt/assertion_levels.hpp"
#include "tdt/graph/compressed-forest.hpp"
#include "tdt/load/forest-compressor.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"

// We store the number of sites without a mutation in af[0], the number of singletons (sites with one mutation) in
// afs[1] and so on. For simplicity we also compute afs[num_samples], that is the number of sites with all samples
// having a derived state.
template <typename AllelicStatePerfectHasher>
class AlleleFrequencySpectrum {
public:
    using value_type     = uint64_t;
    using iterator       = std::vector<value_type>::iterator;
    using const_iterator = std::vector<value_type>::const_iterator;

    AlleleFrequencySpectrum(GenomicSequence& sequence, CompressedForest& forest) {
        auto const num_samples = forest.num_leaves();
        _num_sites             = sequence.num_sites();

        // TODO where to put this line? How do we benchmark it?
        forest.compute_num_samples_below();

        // There are at most \c num_samples many derived samples per site.
        // +1 because there might also be /no/ derived samples.
        _afs.resize(num_samples + 1, 0);

        // For each site, compute how many of the samples have which derived or ancestral state and build the histogram
        // of the number of derived states per site.
        std::array<size_t, AllelicStatePerfectHasher::num_states> num_samples_in_state;
        for (SiteId site = 0; site < _num_sites; ++site) {
            KASSERT(site < sequence.num_sites(), "Site out of bounds", tdt::assert::light);
            num_samples_in_state.fill(0);

            AllelicState const ancestral_state             = sequence.ancestral_state(site);
            auto const         state_of_ancestral_state    = AllelicStatePerfectHasher::to_idx(ancestral_state);
            num_samples_in_state[state_of_ancestral_state] = num_samples;

            for (auto const& mutation: sequence.mutations_at_site(site)) {
                auto const num_samples_below_this_mutation = forest.num_samples_below(mutation.subtree_id());
                auto const state_of_this_mutation = AllelicStatePerfectHasher::to_idx(mutation.allelic_state());

                size_t state_of_parent_mutation = 0;
                // TODO this should be encapsulated in the GenomicSequenceStorage
                if (mutation.parent_mutation_id() == TSK_NULL) {
                    state_of_parent_mutation = state_of_ancestral_state;
                } else {
                    auto parent_mutation =
                        sequence.mutation_by_id(asserting_cast<size_t>(mutation.parent_mutation_id()));
                    state_of_parent_mutation = AllelicStatePerfectHasher::to_idx(parent_mutation.allelic_state());
                }

                KASSERT(
                    num_samples_in_state[state_of_parent_mutation] >= num_samples_below_this_mutation,
                    "There should never be more derived samples than total samples.",
                    tdt::assert::light
                );
                KASSERT(
                    num_samples_in_state[state_of_this_mutation] + num_samples_below_this_mutation <= num_samples,
                    "There should never be more derived samples than total samples.",
                    tdt::assert::light
                );
                num_samples_in_state[state_of_this_mutation] += num_samples_below_this_mutation;
                num_samples_in_state[state_of_parent_mutation] -= num_samples_below_this_mutation;
            }
            for (size_t idx = 0; idx < num_samples_in_state.size(); ++idx) {
                if (idx == state_of_ancestral_state) {
                    continue;
                }
                auto const num_derived_samples = num_samples_in_state[idx];
                KASSERT(
                    num_derived_samples < _afs.size(),
                    "AFS histogram does not have enough bins.",
                    tdt::assert::light
                );
                // TODO Decide on and document what we're doing with sites with no mutations
                if (num_derived_samples != 0) {
                    _afs[num_derived_samples] += 1;
                }
            }
        }
    }

    [[nodiscard]] size_t num_samples() const {
        return _afs.size() - 1;
    }

    [[nodiscard]] size_t num_sites() const {
        return _num_sites;
    }

    [[nodiscard]] value_type frequency(size_t num_derived_state) const {
        return _afs[num_derived_state];
    }

    // For comparison with tskit's results
    [[nodiscard]] auto underlying() const {
        return _afs.data();
    }

    iterator begin() noexcept {
        return _afs.begin();
    }

    iterator end() noexcept {
        return _afs.end();
    }

    const_iterator cbegin() const noexcept {
        return _afs.begin();
    }

    const_iterator cend() const noexcept {
        return _afs.end();
    }

private:
    // TODO Use a different datatype
    std::vector<value_type> _afs;
    size_t                  _num_sites;
};

// inline std::ostream& operator<<(std::ostream& os, AlleleFrequencySpectrum const& afs) {
//     os << "{";
//     for (auto freq_it = afs.cbegin(); freq_it != afs.cend(); ++freq_it) {
//         os << fmt::format("{}, ", *freq_it);
//     }
//     os << "}";
//     return os;
// }
