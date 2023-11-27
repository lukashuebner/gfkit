#pragma once

#include <cstddef>

#include "sfkit/samples/NumSamplesBelowAccessor.hpp"
#include "sfkit/samples/primitives.hpp"
#include "sfkit/sequence/Sequence.hpp"
#include "sfkit/stats/Divergence.hpp"
#include "sfkit/stats/Diversity.hpp"

namespace sfkit::stats {

using sfkit::samples::SampleId;
using sfkit::sequence::SiteId;

class Fst {
public:
    // This is per sequence length, the other statistics are not
    template <typename AlleleFrequencies>
    [[nodiscard]] static double
    fst(SiteId const seq_len, AlleleFrequencies const& allele_freqs_0, AlleleFrequencies const& allele_freqs_1) {
        // For sample sets X and Y, if d(X, Y) is the divergence between X and Y, and d(X) is the diversity of
        // X, then what is computed is $F_{ST} = 1 - 2 * (d(X) + d(Y)) / (d(X) + 2 * d(X, Y) + d(Y))$

        auto const n_0 = allele_freqs_0.num_samples_in_sample_set();
        auto const n_1 = allele_freqs_1.num_samples_in_sample_set();

        auto const d_x  = Diversity::diversity(n_1, allele_freqs_0) / seq_len;
        auto const d_y  = Diversity::diversity(n_1, allele_freqs_1) / seq_len;
        auto const d_xy = Divergence::divergence(n_0, allele_freqs_0, n_1, allele_freqs_1) / seq_len;
        auto const fst  = 1.0 - 2.0 * (d_x + d_y) / (d_x + 2.0 * d_xy + d_y);
        return fst;
    }
};
} // namespace sfkit::stats
