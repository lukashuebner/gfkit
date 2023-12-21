#pragma once

#include "sfkit/samples/primitives.hpp"

namespace sfkit::stats {

using sfkit::samples::SampleId;

class Diversity {
public:
    // TODO Define an AlleleFrequencies concept
    template <typename AlleleFrequencies>
    [[nodiscard]] static double diversity(SampleId num_samples, AlleleFrequencies const& allele_freqs) {
        auto const n  = num_samples;
        double     pi = 0.0;

        allele_freqs.visit(
            // Biallelic visitor
            [&pi, n](auto&& state) {
                auto const n_anc = state.num_ancestral();
                auto const n_der = n - n_anc;
                pi += 2 * static_cast<double>(n_anc * n_der);
            },
            // Multiallelic visitor
            [&pi, n](auto&& state) {
                for (auto const n_state: state) {
                    auto const n_not_state = n - n_state;
                    pi += static_cast<double>(n_state * n_not_state);
                }
            }
        );

        return pi / static_cast<double>(n * (n - 1));
    }
};
} // namespace sfkit::stats
