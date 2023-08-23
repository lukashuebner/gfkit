#pragma once

#include <string>
#include <vector>

#include <kassert/kassert.hpp>
#include <tskit.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/graph/compressed-forest.hpp"
#include "tdt/load/forest-compressor.hpp"
#include "tdt/sequence/allele-frequencies.hpp"
#include "tdt/sequence/allele-frequency-spectrum.hpp"
#include "tdt/sequence/genomic-sequence-storage-factory.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"
#include "tdt/tskit.hpp"
#include "tdt/utils/always_false_v.hpp"

template <typename AllelicStatePerfectHasher = PerfectDNAHasher>
class SuccinctForest {
public:
    using MultiallelicFrequency = typename AlleleFrequencies<AllelicStatePerfectHasher>::MultiallelicFrequency;
    using BiallelicFrequency    = typename AlleleFrequencies<AllelicStatePerfectHasher>::BiallelicFrequency;

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

    AlleleFrequencies<AllelicStatePerfectHasher> allele_frequencies(SampleSet const& sample_set) {
        return AlleleFrequencies<AllelicStatePerfectHasher>(_forest, _sequence, sample_set);
    }

    // TODO Make this const
    [[nodiscard]] double diversity() {
        return diversity(_forest.all_samples());
    }

    [[nodiscard]] double
    diversity(SampleId num_samples, AlleleFrequencies<AllelicStatePerfectHasher> allele_frequencies) {
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

    [[nodiscard]] double divergence(
        SampleId                                     num_samples_1,
        AlleleFrequencies<AllelicStatePerfectHasher> allele_frequencies_1,
        SampleId                                     num_samples_2,
        AlleleFrequencies<AllelicStatePerfectHasher> allele_frequencies_2
    ) {
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

    [[nodiscard]] double divergence(SampleSet const& sample_set_1, SampleSet const& sample_set_2) {
        // TODO Possibly compute both allele frequencies in parallel (high and low order 32 bit)
        auto       allele_freq_1 = allele_frequencies(sample_set_1);
        auto       allele_freq_2 = allele_frequencies(sample_set_2);
        auto const num_samples_1 = sample_set_1.popcount();
        auto const num_samples_2 = sample_set_2.popcount();
        return divergence(num_samples_1, allele_freq_1, num_samples_2, allele_freq_2);
    }

    [[nodiscard]] SiteId num_segregating_sites() {
        return num_segregating_sites(_forest.all_samples());
    }

    // TODO Make this const
    [[nodiscard]] SiteId
    num_segregating_sites(SampleId num_samples, AlleleFrequencies<AllelicStatePerfectHasher> allele_frequencies) {
        SiteId num_segregating_sites = 0;
        allele_frequencies.visit(
            [&num_segregating_sites, num_samples](auto&& state) {
                if (state.num_ancestral() != 0 && state.num_ancestral() != num_samples) {
                    num_segregating_sites++;
                }
            },
            [&num_segregating_sites, num_samples](auto&& state) {
                for (auto const num_samples_in_state: state) {
                    if (num_samples_in_state == num_samples) {
                        return;
                    }
                }
                num_segregating_sites++;
            }
        );

        return num_segregating_sites;
    }

    [[nodiscard]] SiteId num_segregating_sites(SampleSet const& sample_set) {
        auto const num_samples = sample_set.popcount();
        auto const freqs       = allele_frequencies(sample_set);
        return num_segregating_sites(num_samples, freqs);
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
    [[nodiscard]] double fst(SampleSet const& sample_set_1, SampleSet const& sample_set_2) {
        // For sample sets X and Y, if d(X, Y) is the divergence between X and Y, and d(X) is the diversity of X, then
        // what is computed is $F_{ST} = 1 - 2 * (d(X) + d(Y)) / (d(X) + 2 * d(X, Y) + d(Y))$

        auto const n_1             = sample_set_1.popcount();
        auto const n_2             = sample_set_2.popcount();
        auto const allele_freqs_1  = allele_frequencies(sample_set_1);
        auto const allele_freqs_2  = allele_frequencies(sample_set_2);
        auto const sequence_length = _sequence.num_sites();

        auto const d_x  = diversity(n_1, allele_freqs_1) / sequence_length;
        auto const d_y  = diversity(n_2, allele_freqs_2) / sequence_length;
        auto const d_xy = divergence(n_1, allele_freqs_1, n_2, allele_freqs_2) / sequence_length;
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

    // TODO Split up into CompressedForestIO and SequenceForestIO
    // TODO Better naming to distinguish between CompressedForest and SequenceForest
    // friend class CompressedForestIO;
};