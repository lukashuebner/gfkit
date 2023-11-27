#pragma once

#include <cstddef>

#include "sfkit/samples/primitives.hpp"
#include "sfkit/sequence/Sequence.hpp"
#include "sfkit/stats/Diversity.hpp"
#include "sfkit/stats/NumSegregatingSites.hpp"

namespace sfkit::stats {

using sfkit::samples::SampleId;
using sfkit::sequence::SiteId;

class TajimasD {
public:
    template <typename AlleleFrequencies>
    [[nodiscard]] static double tajimas_d(SampleId num_samples, AlleleFrequencies allele_freqs) {
        SampleId const n = num_samples;

        // TODO Does the compiler optimize the following two loops into one?
        double const T = Diversity::diversity(n, allele_freqs);
        double const S = static_cast<double>(NumSegregatingSites::num_segregating_sites(n, allele_freqs));

        double h = 0;
        double g = 0;
        // TODO are there formulas for computing these values more efficiently?
        for (SampleId i = 1; i < num_samples; i++) {
            h += 1.0 / static_cast<double>(i);
            g += 1.0 / static_cast<double>(i * i);
        }
        double const a = (n + 1) / (3 * (n - 1) * h) - 1 / (h * h);
        double const b = 2. * (n * n + n + 3) / (9 * n * (n - 1)) - (n + 2) / (h * n) + g / (h * h);
        double const D = (T - S / h) / sqrt(a * S + (b / (h * h + g)) * S * (S - 1));

        return D;
    }
};
} // namespace sfkit::stats
