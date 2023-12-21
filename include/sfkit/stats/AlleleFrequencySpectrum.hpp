#pragma once

#include <cstddef>
#include <vector>

#include <kassert/kassert.hpp>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/sequence/AlleleFrequencies.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"

namespace sfkit::stats {

using sfkit::samples::SampleId;
using sfkit::sequence::SiteId;
using sfkit::utils::asserting_cast;

// We store the number of sites without a mutation in afs[0], the number of singletons (sites with one mutation) in
// afs[1] and so on. For simplicity we also compute afs[num_samples], that is the number of sites with all samples
// having a derived state.
template <typename AlleleFrequencies>
class AlleleFrequencySpectrum {
public:
    using value_type     = SiteId;
    using iterator       = std::vector<value_type>::iterator;
    using const_iterator = std::vector<value_type>::const_iterator;

    explicit AlleleFrequencySpectrum(AlleleFrequencies const& allele_frequencies) {
        using MultiallelicFrequency = typename AlleleFrequencies::MultiallelicFrequencyT;
        // There are at most \c num_samples many derived samples per site.
        // +1 because there might also be /no/ derived samples.
        auto const num_samples = allele_frequencies.num_samples_in_sample_set();
        _afs.resize(num_samples + 1, 0);

        // For each site, compute how many of the samples have which derived or ancestral state and build the histogram
        // of the number of derived states per site.
        allele_frequencies.visit(
            // Biallelic visitor
            [this, num_samples](auto&& state) {
                auto const num_derived_samples = num_samples - state.num_ancestral();
                KASSERT(
                    num_derived_samples < _afs.size(),
                    "AFS histogram does not have enough bins.",
                    sfkit::assert::light
                );
                // TODO Decide on and document what we're doing with sites with no mutations
                if (num_derived_samples != 0) {
                    _afs[num_derived_samples] += 1;
                }
            },
            // Multiallelic visitor
            [this](auto&& state) {
                for (typename MultiallelicFrequency::Idx state_idx = 0; state_idx < MultiallelicFrequency::num_states;
                     state_idx++) {
                    auto const num_derived_samples = state[state_idx];
                    KASSERT(
                        num_derived_samples < _afs.size(),
                        "AFS histogram does not have enough bins.",
                        sfkit::assert::light
                    );
                    // TODO Decide on and document what we're doing with sites with no mutations
                    if (num_derived_samples != 0 && state_idx != state.ancestral_state_idx()) {
                        _afs[num_derived_samples] += 1;
                    }
                }
            }
        );
    }

    [[nodiscard]] SampleId num_samples() const {
        return asserting_cast<SampleId>(_afs.size()) - 1;
    }

    [[nodiscard]] value_type frequency(size_t num_derived_state) const {
        return _afs[num_derived_state];
    }

    [[nodiscard]] value_type operator[](size_t num_derived_state) const {
        return frequency(num_derived_state);
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

    const_iterator begin() const noexcept {
        return _afs.begin();
    }

    const_iterator end() const noexcept {
        return _afs.end();
    }

private:
    std::vector<value_type> _afs;
};
} // namespace sfkit::stats
