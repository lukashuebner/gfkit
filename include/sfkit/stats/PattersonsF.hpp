#pragma once

#include <variant>

#include "sfkit/samples/NumSamplesBelowAccessor.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/samples/primitives.hpp"

namespace sfkit::stats {

using sfkit::samples::SampleId;

class PattersonsF {
public:
    template <typename AlleleFrequencies>
    [[nodiscard]] static double f2(AlleleFrequencies const& allele_freqs_0, AlleleFrequencies const& allele_freqs_1) {
        using BiallelicFrequency    = typename AlleleFrequencies::BiallelicFrequencyT;
        using MultiallelicFrequency = typename AlleleFrequencies::MultiallelicFrequencyT;

        double const num_samples_0 = allele_freqs_0.num_samples_in_sample_set();
        double const num_samples_1 = allele_freqs_1.num_samples_in_sample_set();
        KASSERT(
            num_samples_0 >= 2.0,
            "We have to draw /two/ samples from the first sample set. It thus must be at least of size 2.",
            sfkit::assert::light
        );
        KASSERT(
            num_samples_1 >= 2.0,
            "We have to draw /two/ samples from the second sample set. It thus must be at least of size 2.",
            sfkit::assert::light
        );

        double f2 = 0.0;

        auto allele_freqs_0_it = allele_freqs_0.cbegin();
        auto allele_freqs_1_it = allele_freqs_1.cbegin();

        while (allele_freqs_0_it != allele_freqs_1.cend()) {
            KASSERT(
                allele_freqs_1_it != allele_freqs_1.cend(),
                "Allele frequency lists have different lengths (different number of sites).",
                sfkit::assert::light
            );

            if (std::holds_alternative<BiallelicFrequency>(*allele_freqs_0_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freqs_1_it)) [[likely]] {
                double const n_anc_0 = std::get<BiallelicFrequency>(*allele_freqs_0_it).num_ancestral();
                double const n_anc_1 = std::get<BiallelicFrequency>(*allele_freqs_1_it).num_ancestral();
                double const n_der_0 = static_cast<double>(num_samples_0) - n_anc_0;
                double const n_der_1 = static_cast<double>(num_samples_1) - n_anc_1;

                f2 += n_anc_0 * (n_anc_0 - 1) * n_der_1 * (n_der_1 - 1) - n_anc_0 * n_der_0 * n_anc_1 * n_der_1;
                f2 += n_der_0 * (n_der_0 - 1) * n_anc_1 * (n_anc_1 - 1) - n_der_0 * n_anc_0 * n_der_1 * n_anc_1;
            } else {
                allele_freqs_0_it.force_multiallelicity();
                allele_freqs_1_it.force_multiallelicity();
                auto const freqs_1_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_0_it);
                auto const freqs_2_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_1_it);

                using Idx = typename MultiallelicFrequency::Idx;
                KASSERT(
                    freqs_1_multiallelic.num_states == freqs_2_multiallelic.num_states,
                    "Allele frequency lists have different lengths (different number of states).",
                    sfkit::assert::light
                );
                for (Idx state = 0; state < freqs_1_multiallelic.num_states; state++) {
                    double const n_state_0     = freqs_1_multiallelic[state];
                    double const n_state_1     = freqs_2_multiallelic[state];
                    double const n_not_state_1 = num_samples_0 - n_state_0;
                    double const n_not_state_2 = num_samples_1 - n_state_1;
                    f2 += n_state_0 * (n_state_0 - 1) * n_not_state_2 * (n_not_state_2 - 1)
                          - n_state_0 * n_not_state_1 * n_state_1 * n_not_state_2;
                }
            }

            allele_freqs_0_it++;
            allele_freqs_1_it++;
        }

        return f2 / (num_samples_0 * (num_samples_0 - 1) * num_samples_1 * (num_samples_1 - 1));
    }

    template <typename AlleleFrequencies>
    [[nodiscard]] static double
    f3(AlleleFrequencies const& allele_freqs_0,
       AlleleFrequencies const& allele_freqs_1,
       AlleleFrequencies const& allele_freqs_2) {
        double f3 = 0.0;

        using MultiallelicFrequency = typename AlleleFrequencies::MultiallelicFrequencyT;
        using BiallelicFrequency    = typename AlleleFrequencies::BiallelicFrequencyT;

        double const num_samples_0 = allele_freqs_0.num_samples_in_sample_set();
        double const num_samples_1 = allele_freqs_1.num_samples_in_sample_set();
        double const num_samples_2 = allele_freqs_2.num_samples_in_sample_set();
        KASSERT(
            num_samples_0 >= 2.0,
            "We have to draw /two/ samples from the first sample set. It thus must be at least of size 2.",
            sfkit::assert::light
        );

        auto allele_freqs_0_it = allele_freqs_0.cbegin();
        auto allele_freqs_1_it = allele_freqs_1.cbegin();
        auto allele_freqs_2_it = allele_freqs_2.cbegin();

        while (allele_freqs_0_it != allele_freqs_0.cend()) {
            KASSERT(
                (allele_freqs_1_it != allele_freqs_1.cend() && allele_freqs_2_it != allele_freqs_2.cend()),
                "Allele frequency lists have different lengths (different number of sites).",
                sfkit::assert::light
            );

            if (std::holds_alternative<BiallelicFrequency>(*allele_freqs_0_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freqs_1_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freqs_2_it)) [[likely]] {
                double const n_anc_0 = std::get<BiallelicFrequency>(*allele_freqs_0_it).num_ancestral();
                double const n_anc_1 = std::get<BiallelicFrequency>(*allele_freqs_1_it).num_ancestral();
                double const n_anc_2 = std::get<BiallelicFrequency>(*allele_freqs_2_it).num_ancestral();
                double const n_der_0 = num_samples_0 - n_anc_0;
                double const n_der_1 = num_samples_1 - n_anc_1;
                double const n_der_2 = num_samples_2 - n_anc_2;

                f3 += n_anc_0 * (n_anc_0 - 1) * n_der_1 * n_der_2 - n_anc_0 * n_der_0 * n_der_1 * n_anc_2;
                f3 += n_der_0 * (n_der_0 - 1) * n_anc_1 * n_anc_2 - n_der_0 * n_anc_0 * n_anc_1 * n_der_2;
            } else {
                allele_freqs_0_it.force_multiallelicity();
                allele_freqs_1_it.force_multiallelicity();
                allele_freqs_2_it.force_multiallelicity();
                auto const freqs_0_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_0_it);
                auto const freqs_1_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_1_it);
                auto const freqs_2_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_2_it);

                using Idx = typename MultiallelicFrequency::Idx;
                KASSERT(
                    (freqs_1_multiallelic.num_states == freqs_1_multiallelic.num_states
                     && freqs_2_multiallelic.num_states == freqs_2_multiallelic.num_states),
                    "Allele frequency lists have different lengths (different number of states).",
                    sfkit::assert::light
                );
                for (Idx state = 0; state < freqs_0_multiallelic.num_states; state++) {
                    double const n_state_0     = freqs_0_multiallelic[state];
                    double const n_state_1     = freqs_1_multiallelic[state];
                    double const n_state_2     = freqs_2_multiallelic[state];
                    double const n_not_state_0 = num_samples_0 - n_state_0;
                    double const n_not_state_1 = num_samples_1 - n_state_1;
                    double const n_not_state_2 = num_samples_2 - n_state_2;
                    f3 += n_state_0 * (n_state_0 - 1) * n_not_state_1 * n_not_state_2
                          - n_state_0 * n_not_state_0 * n_not_state_1 * n_state_2;
                }
            }

            allele_freqs_0_it++;
            allele_freqs_1_it++;
            allele_freqs_2_it++;
        }

        return f3 / (num_samples_0 * (num_samples_0 - 1.0) * num_samples_1 * num_samples_2);
    }

    template <typename AlleleFrequencies>
    [[nodiscard]] static double
    f4(AlleleFrequencies const& allele_freqs_0,
       AlleleFrequencies const& allele_freqs_1,
       AlleleFrequencies const& allele_freqs_2,
       AlleleFrequencies const& allele_freqs_3) {
        double f4 = 0.0;

        using MultiallelicFrequency = typename AlleleFrequencies::MultiallelicFrequencyT;
        using BiallelicFrequency    = typename AlleleFrequencies::BiallelicFrequencyT;

        double const num_samples_0 = allele_freqs_0.num_samples_in_sample_set();
        double const num_samples_1 = allele_freqs_1.num_samples_in_sample_set();
        double const num_samples_2 = allele_freqs_2.num_samples_in_sample_set();
        double const num_samples_3 = allele_freqs_3.num_samples_in_sample_set();

        auto allele_freqs_0_it = allele_freqs_0.cbegin();
        auto allele_freqs_1_it = allele_freqs_1.cbegin();
        auto allele_freqs_2_it = allele_freqs_2.cbegin();
        auto allele_freqs_3_it = allele_freqs_3.cbegin();

        while (allele_freqs_0_it != allele_freqs_1.cend()) {
            KASSERT(
                (allele_freqs_1_it != allele_freqs_1.cend() && allele_freqs_2_it != allele_freqs_2.cend()
                 && allele_freqs_3_it != allele_freqs_3.cend()),
                "Allele frequency lists have different lengths (different number of sites).",
                sfkit::assert::light
            );

            if (std::holds_alternative<BiallelicFrequency>(*allele_freqs_0_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freqs_1_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freqs_2_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freqs_3_it)) [[likely]] {
                double const n_anc_0 = std::get<BiallelicFrequency>(*allele_freqs_0_it).num_ancestral();
                double const n_anc_1 = std::get<BiallelicFrequency>(*allele_freqs_1_it).num_ancestral();
                double const n_anc_2 = std::get<BiallelicFrequency>(*allele_freqs_2_it).num_ancestral();
                double const n_anc_3 = std::get<BiallelicFrequency>(*allele_freqs_3_it).num_ancestral();
                double const n_der_0 = num_samples_0 - n_anc_0;
                double const n_der_1 = num_samples_1 - n_anc_1;
                double const n_der_2 = num_samples_2 - n_anc_2;
                double const n_der_3 = num_samples_3 - n_anc_3;

                // TODO Can we save some computations by rearranging these formulas or does the compiler already
                // do this for us?
                f4 += n_anc_0 * n_der_1 * n_anc_2 * n_der_3 - n_der_0 * n_anc_1 * n_anc_2 * n_der_3;
                f4 += n_der_0 * n_anc_1 * n_der_2 * n_anc_3 - n_anc_0 * n_der_1 * n_der_2 * n_anc_3;
            } else {
                allele_freqs_0_it.force_multiallelicity();
                allele_freqs_1_it.force_multiallelicity();
                allele_freqs_2_it.force_multiallelicity();
                allele_freqs_3_it.force_multiallelicity();
                auto const freqs_0_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_0_it);
                auto const freqs_1_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_1_it);
                auto const freqs_2_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_2_it);
                auto const freqs_3_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_3_it);

                using Idx = typename MultiallelicFrequency::Idx;
                KASSERT(
                    freqs_0_multiallelic.num_states == freqs_1_multiallelic.num_states
                        && freqs_1_multiallelic.num_states == freqs_2_multiallelic.num_states
                        && freqs_2_multiallelic.num_states == freqs_3_multiallelic.num_states,
                    "Allele frequency lists have different lengths (different number of states).",
                    sfkit::assert::light
                );
                for (Idx state = 0; state < freqs_0_multiallelic.num_states; state++) {
                    double const n_state_0 = freqs_0_multiallelic[state];
                    double const n_state_1 = freqs_1_multiallelic[state];
                    double const n_state_2 = freqs_2_multiallelic[state];
                    double const n_state_3 = freqs_3_multiallelic[state];
                    // double const n_not_state_1 = num_samples_1 - n_state_1;
                    double const n_not_state_1 = static_cast<double>(num_samples_1 - n_state_1);
                    double const n_not_state_2 = static_cast<double>(num_samples_2 - n_state_2);
                    double const n_not_state_3 = static_cast<double>(num_samples_3 - n_state_3);
                    f4 += n_state_0 * n_not_state_1 * n_state_2 * n_not_state_3
                          - n_state_0 * n_not_state_1 * n_not_state_2 * n_state_3;
                }
            }

            allele_freqs_0_it++;
            allele_freqs_1_it++;
            allele_freqs_2_it++;
            allele_freqs_3_it++;
        }

        return f4 / (num_samples_0 * num_samples_1 * num_samples_2 * num_samples_3);
    }
};
} // namespace sfkit::stats
