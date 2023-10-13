#pragma once

#include <string>
#include <vector>

#include <kassert/kassert.hpp>
#include <tskit.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/graph/compressed-forest.hpp"
#include "tdt/load/forest-compressor.hpp"
#include "tdt/samples/num-samples-below.hpp"
#include "tdt/sequence/allele-frequencies.hpp"
#include "tdt/sequence/allele-frequency-spectrum.hpp"
#include "tdt/sequence/genomic-sequence-storage-factory.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"
#include "tdt/tskit.hpp"
#include "tdt/utils/always_false_v.hpp"
#include "tdt/utils/tuple_transform.hpp"

template <typename AllelicStatePerfectHasher = PerfectDNAHasher>
class SuccinctForest {
public:
    SuccinctForest(
        TSKitTreeSequence&& tree_sequence, CompressedForest&& compressed_forest, GenomicSequence&& genomic_sequence
    )
        : _tree_sequence(std::move(tree_sequence)),
          _forest(compressed_forest),
          _sequence(genomic_sequence) {}

    SuccinctForest(tsk_treeseq_t&& ts_tree_sequence) : SuccinctForest(TSKitTreeSequence(ts_tree_sequence)) {}

    SuccinctForest(TSKitTreeSequence&& tree_sequence) : _tree_sequence(std::move(tree_sequence)) {
        ForestCompressor       forest_compressor(_tree_sequence);
        GenomicSequenceFactory sequence_factory(_tree_sequence);
        _forest   = forest_compressor.compress(sequence_factory);
        _sequence = sequence_factory.move_storage();
    }

    auto allele_frequencies(SampleSet const& samples) {
        // return AlleleFrequencies<AllelicStatePerfectHasher>(_forest, _sequence, samples);
        auto num_samples_below = NumSamplesBelowAccessor(
            std::make_shared<NumSamplesBelow<1>>(
                _forest.postorder_edges(),
                std::array<std::reference_wrapper<SampleSet const>, 1>{std::cref(samples)}
            ),
            0
        );
        return allele_frequencies(num_samples_below);
    }

    template <NumSamplesBelowAccessorC NumSamplesBelowAccessorT>
    auto allele_frequencies(NumSamplesBelowAccessorT const& num_samples_below) {
        return AlleleFrequencies<AllelicStatePerfectHasher, NumSamplesBelowAccessorT>(
            _forest,
            _sequence,
            num_samples_below
        );
    }

    auto allele_frequencies(SampleSet const& samples_0, SampleSet const& samples_1) {
        return tuple_transform(
            [this](auto&& e) { return allele_frequencies(e); },
            NumSamplesBelowFactory::build(_forest.postorder_edges(), samples_0, samples_1)
        );
    }

    auto allele_frequencies(SampleSet const& samples_0, SampleSet const& samples_1, SampleSet const& samples_2) {
        SampleSet empty_sample_set = SampleSet(samples_0.overall_num_samples());
        auto [num_samples_below_0, num_samples_below_1, num_samples_below_2, dummy] =
            NumSamplesBelowFactory::build(_forest.postorder_edges(), samples_0, samples_1, samples_2, empty_sample_set);
        return std::tuple(
            allele_frequencies(num_samples_below_0),
            allele_frequencies(num_samples_below_1),
            allele_frequencies(num_samples_below_2)
        );
    }

    auto allele_frequencies(
        SampleSet const& samples_0, SampleSet const& samples_1, SampleSet const& samples_2, SampleSet const& samples_3
    ) {
        return tuple_transform(
            [this](auto&& e) { return allele_frequencies(e); },
            NumSamplesBelowFactory::build(_forest.postorder_edges(), samples_0, samples_1, samples_2, samples_3)
        );
    }

    // TODO Make this const
    [[nodiscard]] double diversity() {
        return diversity(_forest.all_samples());
    }

    template <
        typename AlleleFrequenciesT =
            AlleleFrequencies<AllelicStatePerfectHasher, NumSamplesBelowAccessor<NumSamplesBelow<1>>>>
    [[nodiscard]] double diversity(SampleId num_samples, AlleleFrequenciesT allele_frequencies) {
        auto const n  = num_samples;
        double     pi = 0.0;

        allele_frequencies.visit(
            // Biallelic visitor
            [&pi, n](auto&& state) {
                auto const n_anc = state.num_ancestral();
                auto const n_der = n - n_anc;
                pi += 2 * static_cast<double>(n_anc * n_der);
            },
            // Multiallelic visitor
            [&pi, n](auto&& state) {
                // TODO Does the compiler unroll this?
                for (auto const n_state: state) {
                    auto const n_not_state = n - n_state;
                    pi += static_cast<double>(n_state * n_not_state);
                }
            }
        );

        return pi / static_cast<double>(n * (n - 1));
    }

    [[nodiscard]] double diversity(SampleSet const& sample_set) {
        SampleId const num_samples = sample_set.popcount();
        auto const     freqs       = allele_frequencies(sample_set);
        return diversity(num_samples, freqs);
    }

    [[nodiscard]] AlleleFrequencySpectrum<AllelicStatePerfectHasher> allele_frequency_spectrum() {
        return AlleleFrequencySpectrum<AllelicStatePerfectHasher>(allele_frequencies(_forest.all_samples()));
    }

    [[nodiscard]] AlleleFrequencySpectrum<AllelicStatePerfectHasher> allele_frequency_spectrum(SampleSet sample_set) {
        return AlleleFrequencySpectrum<AllelicStatePerfectHasher>(allele_frequencies(sample_set));
    }

    template <
        typename AlleleFrequenciesT =
            AlleleFrequencies<AllelicStatePerfectHasher, NumSamplesBelowAccessor<NumSamplesBelow<1>>>>
    [[nodiscard]] double divergence(
        SampleId           num_samples_1,
        AlleleFrequenciesT allele_frequencies_1,
        SampleId           num_samples_2,
        AlleleFrequenciesT allele_frequencies_2
    ) {
        using MultiallelicFrequency = typename AlleleFrequenciesT::MultiallelicFrequency;
        using BiallelicFrequency    = typename AlleleFrequenciesT::BiallelicFrequency;

        double divergence       = 0.0;
        auto   allele_freq_1_it = allele_frequencies_1.cbegin();
        auto   allele_freq_2_it = allele_frequencies_2.cbegin();
        while (allele_freq_1_it != allele_frequencies_1.cend()) {
            KASSERT(
                allele_freq_2_it != allele_frequencies_2.cend(),
                "Allele frequency lists have different lengths (different number of sites).",
                tdt::assert::light
            );

            if (std::holds_alternative<BiallelicFrequency>(*allele_freq_1_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freq_2_it)) [[likely]] {
                double const n_anc_1 = std::get<BiallelicFrequency>(*allele_freq_1_it).num_ancestral();
                double const n_anc_2 = std::get<BiallelicFrequency>(*allele_freq_2_it).num_ancestral();
                double const n_der_1 = num_samples_1 - n_anc_1;
                double const n_der_2 = num_samples_2 - n_anc_2;
                divergence += static_cast<double>(n_anc_1 * n_der_2 + n_der_1 * n_anc_2);
            } else {
                allele_freq_1_it.force_multiallelicity();
                allele_freq_2_it.force_multiallelicity();
                auto const freq_1_multiallelic = std::get<MultiallelicFrequency>(*allele_freq_1_it);
                auto const freq_2_multiallelic = std::get<MultiallelicFrequency>(*allele_freq_2_it);

                using Idx = typename MultiallelicFrequency::Idx;
                // TODO Check if the compiler unrolls this
                for (Idx state = 0; state < freq_1_multiallelic.num_states; state++) {
                    // The number of samples in sample set 2 which are not in
                    double const n_state_1     = freq_1_multiallelic[state];
                    double const n_not_state_2 = num_samples_2 - freq_2_multiallelic[state];
                    divergence += n_state_1 * n_not_state_2;
                }
            }

            allele_freq_1_it++;
            allele_freq_2_it++;
        }

        return divergence / (static_cast<double>(num_samples_1 * num_samples_2));
    }

    [[nodiscard]] double divergence(SampleSet const& sample_set_0, SampleSet const& sample_set_1) {
        auto [allele_freqs_0, allele_freqs_1] = allele_frequencies(sample_set_0, sample_set_1);
        auto const num_samples_0              = sample_set_0.popcount();
        auto const num_samples_1              = sample_set_1.popcount();
        return divergence(num_samples_0, allele_freqs_0, num_samples_1, allele_freqs_1);
    }

    [[nodiscard]] double f2(SampleSet const& sample_set_0, SampleSet const& sample_set_1) {
        auto [allele_freqs_0, allele_freqs_1] = allele_frequencies(sample_set_0, sample_set_1);
        using BiallelicFrequency              = typename decltype(allele_freqs_0)::BiallelicFrequency;
        using MultiallelicFrequency           = typename decltype(allele_freqs_0)::MultiallelicFrequency;

        double const num_samples_0 = sample_set_0.popcount();
        double const num_samples_1 = sample_set_1.popcount();
        KASSERT(
            num_samples_0 >= 2.0,
            "We have to draw /two/ samples from the first sample set. It thus must be at least of size 2.",
            tdt::assert::light
        );
        KASSERT(
            num_samples_1 >= 2.0,
            "We have to draw /two/ samples from the second sample set. It thus must be at least of size 2.",
            tdt::assert::light
        );

        double f2 = 0.0;

        auto allele_freqs_0_it = allele_freqs_0.cbegin();
        auto allele_freqs_1_it = allele_freqs_1.cbegin();

        while (allele_freqs_0_it != allele_freqs_1.cend()) {
            KASSERT(
                allele_freqs_1_it != allele_freqs_1.cend(),
                "Allele frequency lists have different lengths (different number of sites).",
                tdt::assert::light
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
                // TODO Check if the compiler unrolls this
                KASSERT(
                    freqs_1_multiallelic.num_states == freqs_2_multiallelic.num_states,
                    "Allele frequency lists have different lengths (different number of states).",
                    tdt::assert::light
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

    [[nodiscard]] double
    f3(SampleSet const& sample_set_0, SampleSet const& sample_set_1, SampleSet const& sample_set_2) {
        double f3 = 0.0;

        auto const [allele_freqs_1, allele_freqs_2, allele_freqs_3] =
            allele_frequencies(sample_set_0, sample_set_1, sample_set_2);

        using MultiallelicFrequency = typename decltype(allele_freqs_1)::MultiallelicFrequency;
        using BiallelicFrequency    = typename decltype(allele_freqs_1)::BiallelicFrequency;

        double const num_samples_1 = sample_set_0.popcount();
        double const num_samples_2 = sample_set_1.popcount();
        double const num_samples_3 = sample_set_2.popcount();
        KASSERT(
            num_samples_1 >= 2.0,
            "We have to draw /two/ samples from the first sample set. It thus must be at least of size 2.",
            tdt::assert::light
        );

        auto allele_freqs_1_it = allele_freqs_1.cbegin();
        auto allele_freqs_2_it = allele_freqs_2.cbegin();
        auto allele_freqs_3_it = allele_freqs_3.cbegin();

        while (allele_freqs_1_it != allele_freqs_1.cend()) {
            KASSERT(
                (allele_freqs_2_it != allele_freqs_2.cend() && allele_freqs_3_it != allele_freqs_3.cend()),
                "Allele frequency lists have different lengths (different number of sites).",
                tdt::assert::light
            );

            if (std::holds_alternative<BiallelicFrequency>(*allele_freqs_1_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freqs_2_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freqs_3_it)) [[likely]] {
                double const n_anc_1 = std::get<BiallelicFrequency>(*allele_freqs_1_it).num_ancestral();
                double const n_anc_2 = std::get<BiallelicFrequency>(*allele_freqs_2_it).num_ancestral();
                double const n_anc_3 = std::get<BiallelicFrequency>(*allele_freqs_3_it).num_ancestral();
                double const n_der_1 = num_samples_1 - n_anc_1;
                double const n_der_2 = num_samples_2 - n_anc_2;
                double const n_der_3 = num_samples_3 - n_anc_3;

                f3 += n_anc_1 * (n_anc_1 - 1) * n_der_2 * n_der_3 - n_anc_1 * n_der_1 * n_der_2 * n_anc_3;
                f3 += n_der_1 * (n_der_1 - 1) * n_anc_2 * n_anc_3 - n_der_1 * n_anc_1 * n_anc_2 * n_der_3;
            } else {
                allele_freqs_1_it.force_multiallelicity();
                allele_freqs_2_it.force_multiallelicity();
                allele_freqs_3_it.force_multiallelicity();
                auto const freqs_1_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_1_it);
                auto const freqs_2_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_2_it);
                auto const freqs_3_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_3_it);

                using Idx = typename MultiallelicFrequency::Idx;
                // TODO Check if the compiler unrolls this
                KASSERT(
                    (freqs_1_multiallelic.num_states == freqs_2_multiallelic.num_states
                     && freqs_2_multiallelic.num_states == freqs_3_multiallelic.num_states),
                    "Allele frequency lists have different lengths (different number of states).",
                    tdt::assert::light
                );
                for (Idx state = 0; state < freqs_1_multiallelic.num_states; state++) {
                    double const n_state_1     = freqs_1_multiallelic[state];
                    double const n_state_2     = freqs_2_multiallelic[state];
                    double const n_state_3     = freqs_3_multiallelic[state];
                    double const n_not_state_1 = num_samples_1 - n_state_1;
                    double const n_not_state_2 = num_samples_2 - n_state_2;
                    double const n_not_state_3 = num_samples_3 - n_state_3;
                    f3 += n_state_1 * (n_state_1 - 1) * n_not_state_2 * n_not_state_3
                          - n_state_1 * n_not_state_1 * n_not_state_2 * n_state_3;
                }
            }

            allele_freqs_1_it++;
            allele_freqs_2_it++;
            allele_freqs_3_it++;
        }

        return f3 / (num_samples_1 * (num_samples_1 - 1.0) * num_samples_2 * num_samples_3);
    }

    [[nodiscard]] double
    f4(SampleSet const& sample_set_1,
       SampleSet const& sample_set_2,
       SampleSet const& sample_set_3,
       SampleSet const& sample_set_4) {
        double f4 = 0.0;

        auto const [allele_freqs_1, allele_freqs_2, allele_freqs_3, allele_freqs_4] =
            allele_frequencies(sample_set_1, sample_set_2, sample_set_3, sample_set_4);
        using MultiallelicFrequency = typename decltype(allele_freqs_1)::MultiallelicFrequency;
        using BiallelicFrequency    = typename decltype(allele_freqs_1)::BiallelicFrequency;

        double const num_samples_1 = sample_set_1.popcount();
        double const num_samples_2 = sample_set_2.popcount();
        double const num_samples_3 = sample_set_3.popcount();
        double const num_samples_4 = sample_set_4.popcount();

        auto allele_freqs_1_it = allele_freqs_1.cbegin();
        auto allele_freqs_2_it = allele_freqs_2.cbegin();
        auto allele_freqs_3_it = allele_freqs_3.cbegin();
        auto allele_freqs_4_it = allele_freqs_4.cbegin();

        while (allele_freqs_1_it != allele_freqs_1.cend()) {
            KASSERT(
                (allele_freqs_2_it != allele_freqs_2.cend() && allele_freqs_3_it != allele_freqs_3.cend()
                 && allele_freqs_4_it != allele_freqs_4.cend()),
                "Allele frequency lists have different lengths (different number of sites).",
                tdt::assert::light
            );

            if (std::holds_alternative<BiallelicFrequency>(*allele_freqs_1_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freqs_2_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freqs_3_it)
                && std::holds_alternative<BiallelicFrequency>(*allele_freqs_4_it)) [[likely]] {
                double const n_anc_1 = std::get<BiallelicFrequency>(*allele_freqs_1_it).num_ancestral();
                double const n_anc_2 = std::get<BiallelicFrequency>(*allele_freqs_2_it).num_ancestral();
                double const n_anc_3 = std::get<BiallelicFrequency>(*allele_freqs_3_it).num_ancestral();
                double const n_anc_4 = std::get<BiallelicFrequency>(*allele_freqs_4_it).num_ancestral();
                double const n_der_1 = num_samples_1 - n_anc_1;
                double const n_der_2 = num_samples_2 - n_anc_2;
                double const n_der_3 = num_samples_3 - n_anc_3;
                double const n_der_4 = num_samples_4 - n_anc_4;

                // TODO Can we save some computations by rearranging these formulas or does the compiler already do this
                // for us?
                f4 += n_anc_1 * n_der_2 * n_anc_3 * n_der_4 - n_der_1 * n_anc_2 * n_anc_3 * n_der_4;
                f4 += n_der_1 * n_anc_2 * n_der_3 * n_anc_4 - n_anc_1 * n_der_2 * n_der_3 * n_anc_4;
            } else {
                allele_freqs_1_it.force_multiallelicity();
                allele_freqs_2_it.force_multiallelicity();
                allele_freqs_3_it.force_multiallelicity();
                allele_freqs_4_it.force_multiallelicity();
                auto const freqs_1_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_1_it);
                auto const freqs_2_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_2_it);
                auto const freqs_3_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_3_it);
                auto const freqs_4_multiallelic = std::get<MultiallelicFrequency>(*allele_freqs_4_it);

                using Idx = typename MultiallelicFrequency::Idx;
                // TODO Check if the compiler unrolls this
                KASSERT(
                    freqs_1_multiallelic.num_states == freqs_2_multiallelic.num_states
                        && freqs_2_multiallelic.num_states == freqs_3_multiallelic.num_states
                        && freqs_3_multiallelic.num_states == freqs_4_multiallelic.num_states,
                    "Allele frequency lists have different lengths (different number of states).",
                    tdt::assert::light
                );
                for (Idx state = 0; state < freqs_1_multiallelic.num_states; state++) {
                    double const n_state_1 = freqs_1_multiallelic[state];
                    double const n_state_2 = freqs_2_multiallelic[state];
                    double const n_state_3 = freqs_3_multiallelic[state];
                    double const n_state_4 = freqs_4_multiallelic[state];
                    // double const n_not_state_1 = num_samples_1 - n_state_1;
                    double const n_not_state_2 = static_cast<double>(num_samples_2 - n_state_2);
                    double const n_not_state_3 = static_cast<double>(num_samples_3 - n_state_3);
                    double const n_not_state_4 = static_cast<double>(num_samples_4 - n_state_4);
                    f4 += n_state_1 * n_not_state_2 * n_state_3 * n_not_state_4
                          - n_state_1 * n_not_state_2 * n_not_state_3 * n_state_4;
                }
            }

            allele_freqs_1_it++;
            allele_freqs_2_it++;
            allele_freqs_3_it++;
            allele_freqs_4_it++;
        }

        return f4 / (num_samples_1 * num_samples_2 * num_samples_3 * num_samples_4);
    }

    // TODO Make this const
    [[nodiscard]] SiteId
    num_segregating_sites(SampleId num_samples, AlleleFrequencies<AllelicStatePerfectHasher> allele_frequencies) {
        // We compute the number of segregating sites in this way in order to be compatible with tskit's
        // definition. See https://doi.org/10.1534/genetics.120.303253 and tskit's trees.c
        size_t num_segregating_sites = 0;
        allele_frequencies.visit(
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
                KASSERT(num_states > 0ul, "There are no allelic states at this site.", tdt::assert::light);
                num_segregating_sites += num_states - 1;
            }
        );

        return asserting_cast<SiteId>(num_segregating_sites);
    }

    [[nodiscard]] SiteId num_segregating_sites(SampleSet const& sample_set) {
        auto const num_samples = sample_set.popcount();
        auto const freqs       = allele_frequencies(sample_set);
        return num_segregating_sites(num_samples, freqs);
    }

    [[nodiscard]] SiteId num_segregating_sites() {
        return num_segregating_sites(_forest.all_samples());
    }

    [[nodiscard]] double tajimas_d() {
        SampleId const n = num_samples();

        auto const allele_freqs = allele_frequencies(_forest.all_samples());
        // TODO Does the compiler optimize the following two loops into one?
        double const T = diversity(n, allele_freqs);
        double const S = static_cast<double>(num_segregating_sites(n, allele_freqs));

        double h = 0;
        double g = 0;
        // TODO are there formulas for computing these values more efficiently?
        for (SampleId i = 1; i < num_samples(); i++) {
            h += 1.0 / static_cast<double>(i);
            g += 1.0 / static_cast<double>(i * i);
        }
        double const a = (n + 1) / (3 * (n - 1) * h) - 1 / (h * h);
        double const b = 2. * (n * n + n + 3) / (9 * n * (n - 1)) - (n + 2) / (h * n) + g / (h * h);
        double const D = (T - S / h) / sqrt(a * S + (b / (h * h + g)) * S * (S - 1));

        return D;
    }

    // This is per sequence length, the other statistics are not
    [[nodiscard]] double fst(SampleSet const& sample_set_0, SampleSet const& sample_set_1) {
        // For sample sets X and Y, if d(X, Y) is the divergence between X and Y, and d(X) is the diversity of X,
        // then what is computed is $F_{ST} = 1 - 2 * (d(X) + d(Y)) / (d(X) + 2 * d(X, Y) + d(Y))$

        auto const n_0                        = sample_set_0.popcount();
        auto const n_1                        = sample_set_1.popcount();
        auto [allele_freqs_0, allele_freqs_1] = allele_frequencies(sample_set_0, sample_set_1);
        auto const sequence_length            = _sequence.num_sites();

        auto const d_x  = diversity(n_1, allele_freqs_0) / sequence_length;
        auto const d_y  = diversity(n_1, allele_freqs_1) / sequence_length;
        auto const d_xy = divergence(n_0, allele_freqs_0, n_1, allele_freqs_1) / sequence_length;
        auto const fst  = 1.0 - 2.0 * (d_x + d_y) / (d_x + 2.0 * d_xy + d_y);
        return fst;
    }

    [[nodiscard]] SiteId num_sites() const {
        return _sequence.num_sites();
    }

    [[nodiscard]] SampleId num_samples() const {
        return _forest.num_samples();
    }

    [[nodiscard]] SampleSet all_samples() const {
        return _forest.all_samples();
    }

    [[nodiscard]] GenomicSequence const& sequence() const {
        return _sequence;
    }

    [[nodiscard]] CompressedForest const& forest() const {
        return _forest;
    }

    [[nodiscard]] TreeId num_trees() const {
        return _forest.num_trees();
    }

    [[nodiscard]] MutationId num_mutations() const {
        return _sequence.num_mutations();
    }

    [[nodiscard]] auto num_unique_subtrees() const {
        return _forest.num_unique_subtrees();
    }

    [[nodiscard]] auto num_subtrees_with_mutations() const {
        return _sequence.subtrees_with_mutations().size();
    }

    [[nodiscard]] TSKitTreeSequence& tree_sequence() {
        return _tree_sequence;
    }

    // TODO We should not store this!!
    [[nodiscard]] TSKitTreeSequence const& tree_sequence() const {
        return _tree_sequence;
    }

private:
    TSKitTreeSequence _tree_sequence;
    CompressedForest  _forest;
    GenomicSequence   _sequence;

    // TODO Split up into CompressedForestIO and SuccinctForestIO
    // TODO Better naming to distinguish between CompressedForest and SuccinctForest
    // friend class CompressedForestIO;
};
