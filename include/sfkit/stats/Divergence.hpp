#pragma once

#include <cstddef>
#include <variant>

#include "sfkit/samples/primitives.hpp"
#include "sfkit/sequence/Sequence.hpp"

namespace sfkit::stats {

using sfkit::samples::SampleId;
using sfkit::sequence::SiteId;

class Divergence {
public:
    template <typename AlleleFrequencies>
    [[nodiscard]] static double divergence(
        SampleId          num_samples_0,
        AlleleFrequencies allele_freqs_0,
        SampleId          num_samples_1,
        AlleleFrequencies allele_freqs_1
    ) {
        using MultiallelicFrequency = typename AlleleFrequencies::MultiallelicFrequency;
        using BiallelicFrequency    = typename AlleleFrequencies::BiallelicFrequency;

        double divergence       = 0.0;
        auto   allele_freq_0_it = allele_freqs_0.cbegin();
        auto   allele_freq_1_it = allele_freqs_1.cbegin();
        while (allele_freq_0_it != allele_freqs_0.cend()) {
            KASSERT(
                allele_freq_1_it != allele_freqs_1.cend(),
                "Allele frequency lists have different lengths (different number of sites).",
                sfkit::assert::light
            );

            if (std::holds_alternative<BiallelicFrequency>(*allele_freq_0_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freq_1_it)) [[likely]] {
                double const n_anc_0 = std::get<BiallelicFrequency>(*allele_freq_0_it).num_ancestral();
                double const n_anc_1 = std::get<BiallelicFrequency>(*allele_freq_1_it).num_ancestral();
                double const n_der_0 = num_samples_0 - n_anc_0;
                double const n_der_1 = num_samples_1 - n_anc_1;
                divergence += static_cast<double>(n_anc_0 * n_der_1 + n_der_0 * n_anc_1);
            } else {
                allele_freq_0_it.force_multiallelicity();
                allele_freq_1_it.force_multiallelicity();
                auto const freq_0_multiallelic = std::get<MultiallelicFrequency>(*allele_freq_0_it);
                auto const freq_1_multiallelic = std::get<MultiallelicFrequency>(*allele_freq_1_it);

                using Idx = typename MultiallelicFrequency::Idx;
                for (Idx state = 0; state < freq_0_multiallelic.num_states; state++) {
                    double const n_state_0     = freq_0_multiallelic[state];
                    double const n_not_state_1 = num_samples_1 - freq_1_multiallelic[state];
                    divergence += n_state_0 * n_not_state_1;
                }
            }

            allele_freq_0_it++;
            allele_freq_1_it++;
        }

        return divergence / (static_cast<double>(num_samples_0 * num_samples_1));
    }
};
} // namespace sfkit::stats
