#pragma once

#include <cstddef>

#include "sfkit/samples/primitives.hpp"
#include "sfkit/sequence/Sequence.hpp"

namespace sfkit::stats {

using sfkit::samples::SampleId;
using sfkit::sequence::SiteId;
using sfkit::utils::asserting_cast;

class NumSegregatingSites {
public:
    template <typename AlleleFrequencies>
    [[nodiscard]] static SiteId num_segregating_sites(SampleId num_samples, AlleleFrequencies const& allele_freqs) {
        // We compute the number of segregating sites in this way in order to be compatible with tskit's
        // definition. See https://doi.org/10.1534/genetics.120.303253 and tskit's trees.c
        size_t num_segregating_sites = 0;
        allele_freqs.visit(
            [&num_segregating_sites, num_samples](auto&& state) {
                num_segregating_sites += (state.num_ancestral() > 0 && state.num_ancestral() < num_samples);
            },
            [&num_segregating_sites, num_samples](auto&& states) {
                size_t num_states = 0;
                for (auto const num_samples_in_state: states) {
                    if (num_samples_in_state > 0ul) {
                        num_states++;
                    }
                }
                KASSERT(num_states > 0ul, "There are no allelic states at this site.", sfkit::assert::light);
                num_segregating_sites += num_states - 1;
            }
        );

        return asserting_cast<SiteId>(num_segregating_sites);
    }
};
} // namespace sfkit::stats
