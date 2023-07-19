#pragma once

#include <string>
#include <vector>

#include <kassert/kassert.hpp>
#include <tskit.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/graph/compressed-forest.hpp"
#include "tdt/load/forest-compressor.hpp"
#include "tdt/sequence/allele-frequencies.hpp"
#include "tdt/sequence/genomic-sequence-storage-factory.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"
#include "tdt/tskit.hpp"

// TODO encapsulate CompresedForest and GenomicSequence functions in this class
class SequenceForest {
public:
    // TODO Rename to GenomicSequenceStore
    SequenceForest(
        TSKitTreeSequence&& tree_sequence, CompressedForest&& compressed_forest, GenomicSequence&& genomic_sequence
    )
        : _tree_sequence(std::move(tree_sequence)),
          _forest(compressed_forest),
          _sequence(genomic_sequence) {}

    SequenceForest(tsk_treeseq_t&& ts_tree_sequence) : SequenceForest(TSKitTreeSequence(ts_tree_sequence)) {}

    SequenceForest(TSKitTreeSequence&& tree_sequence) : _tree_sequence(std::move(tree_sequence)) {
        ForestCompressor       forest_compressor(_tree_sequence);
        GenomicSequenceFactory sequence_factory(_tree_sequence);
        _forest   = forest_compressor.compress(sequence_factory);
        _sequence = sequence_factory.move_storage();
    }

    AlleleFrequencies allele_frequencies(SampleSet const& sample_set) {
        return AlleleFrequencies(_forest, _sequence, sample_set);
    }

    // TODO Make this const
    [[nodiscard]] double diversity() {
        return diversity(_forest.all_samples());
    }

    [[nodiscard]] double diversity(SampleSet const& sample_set) {
        auto const n  = num_samples();
        double     pi = 0.0;
        for (uint64_t frequency: allele_frequencies(sample_set)) {
            pi += static_cast<double>(frequency * (n - frequency));
        }

        return pi / (0.5 * static_cast<double>(n * (n - 1)));
    }

    [[nodiscard]] double divergence(SampleSet const& sample_set_1, SampleSet const& sample_set_2) {
        auto const n1 = sample_set_1.popcount();
        auto const n2 = sample_set_2.popcount();

        double divergence = 0.0;
        // TODO Combine both frequencies computations into one
        auto allele_freq_1    = allele_frequencies(sample_set_1);
        auto allele_freq_2    = allele_frequencies(sample_set_2);
        auto allele_freq_1_it = allele_freq_1.cbegin();
        auto allele_freq_2_it = allele_freq_2.cbegin();
        while (allele_freq_1_it != allele_freq_1.cend()) {
            KASSERT(
                allele_freq_2_it != allele_freq_2.cend(),
                "Allele frequency lists have different lengths (different number of sites).",
                tdt::assert::light
            );
            auto const freq1 = *allele_freq_1_it;
            auto const freq2 = *allele_freq_2_it;
            // divergence += static_cast<double>(freq1 * (n2 - freq2));
            divergence += static_cast<double>(freq1 * (n2 - freq2) + freq2 * (n1 - freq1));
            allele_freq_1_it++;
            allele_freq_2_it++;
        }

        return divergence / (static_cast<double>(n1 * n2));
    }

    [[nodiscard]] uint64_t num_segregating_sites() {
        return num_segregating_sites(_forest.all_samples());
    }

    // TODO Make this const
    [[nodiscard]] uint64_t num_segregating_sites(SampleSet const& sample_set) {
        size_t num_segregating_sites = 0;
        // TODO Pass sample set as shared_ptr? 
        for (uint64_t frequency: allele_frequencies(sample_set)) {
            if (frequency != 0 && frequency != num_samples()) {
                num_segregating_sites++;
            }
        }
        return num_segregating_sites;
    }

    [[nodiscard]] double tajimas_d() {
        double const n = static_cast<double>(num_samples());
        // TODO Does the compiler optimize the following two loops into one?
        double const T = diversity();
        double const S = static_cast<double>(num_segregating_sites());

        double h = 0;
        double g = 0;
        // TODO are there formulas for computing these values more efficiently?
        for (uint64_t i = 1; i < num_samples(); i++) {
            h += 1.0 / static_cast<double>(i);
            g += 1.0 / static_cast<double>(i * i);
        }
        double const a = (n + 1) / (3 * (n - 1) * h) - 1 / (h * h);
        double const b = 2 * (n * n + n + 3) / (9 * n * (n - 1)) - (n + 2) / (h * n) + g / (h * h);
        double const D = (T - S / h) / sqrt(a * S + (b / (h * h + g)) * S * (S - 1));

        return D;
    }

    [[nodiscard]] uint64_t num_sites() const {
        return _sequence.num_sites();
    }

    [[nodiscard]] uint64_t num_samples() const {
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

    // TODO We should not store this!!
    [[nodiscard]] TSKitTreeSequence const& tree_sequence() const {
        return _tree_sequence;
    }

private:
    TSKitTreeSequence _tree_sequence;
    CompressedForest  _forest;
    GenomicSequence   _sequence;
};
