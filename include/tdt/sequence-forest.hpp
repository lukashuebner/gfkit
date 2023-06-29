#pragma once

#include <string>

#include <kassert/kassert.hpp>
#include <tskit.h>

#include "assertion_levels.hpp"
#include "graph/compressed-forest.hpp"
#include "load/forest-compressor.hpp"
#include "sequence/genomic-sequence-storage-factory.hpp"
#include "sequence/genomic-sequence-storage.hpp"
#include "tskit.hpp"

// TODO This should reside in its own file
// TODO This class shares functionality with the AlleleFrequencySpectrum class. We should refactor this.
// The historical reason for this duplication is that AlleleFrequencySpectrum was written under the assumption of
// multiallelicity and this class was written under the assumption of biallelicity when I implemented the diversity
// statistics.
class AlleleFrequencies {
public:
    AlleleFrequencies(CompressedForest& compressed_forest, GenomicSequenceStorage const& sequence_store)
        : _forest(compressed_forest),
          _sequence(sequence_store) {}

    auto begin() const {
        return allele_frequency_iterator{_forest, _sequence};
    }

    auto end() const {
        return allele_frequency_iterator::sentinel{};
    }

    auto cbegin() const {
        return allele_frequency_iterator{_forest, _sequence};
    }

    auto cend() const {
        return allele_frequency_iterator::sentinel{};
    }

    class allele_frequency_iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = uint64_t;
        using pointer           = value_type*;
        using reference         = value_type&;
        struct sentinel {};

        allele_frequency_iterator(CompressedForest& compressed_forest, GenomicSequenceStorage const& sequence_store)
            : _forest(compressed_forest),
              _sequence(sequence_store),
              _site(0) {
            // TODO where to put this line? Cache the sample computation?
            compressed_forest.compute_num_samples_below();
            _update_state();
        }

        auto num_samples() const {
            return _forest.num_samples();
        }

        auto num_sites() const {
            return _sequence.num_sites();
        }

        allele_frequency_iterator& operator++() {
            _site++;
            if (_site < num_sites()) {
                _update_state();
            }
            return *this;
        }

        allele_frequency_iterator operator++(int) {
            allele_frequency_iterator tmp(*this);
            this->                    operator++();
            return tmp;
        }

        bool operator==(allele_frequency_iterator const& other) const {
            return &_forest == &other._forest && &_sequence == &other._sequence;
        }

        bool operator!=(allele_frequency_iterator const& other) const {
            return !(*this == other);
        }

        bool operator==(sentinel) {
            return _site == num_sites();
        }

        reference operator*() {
            return _frequency;
        }

        pointer operator->() {
            return &_frequency;
        }

    private:
        CompressedForest&             _forest;
        GenomicSequenceStorage const& _sequence;
        SiteId                        _site;
        uint64_t                      _frequency;

        void _update_state() {
            // TODO Generalize this to multiallelic sites
            KASSERT(_site < num_sites(), "Site index out of bounds", tdt::assert::light);

            uint64_t           num_samples_in_ancestral_state = num_samples();
            AllelicState const ancestral_state                = _sequence.ancestral_state(_site);

#if KASSERT_ENABLED(TDT_ASSERTION_LEVEL_NORMAL)
            AllelicState derived_state = InvalidAllelicState;
#endif
            for (auto const& mutation: _sequence.mutations_at_site(_site)) {
                auto const num_samples_below_this_mutation = _forest.num_samples_below(mutation.subtree_id());
                auto const this_mutations_state            = mutation.allelic_state();

#if KASSERT_ENABLED(TDT_ASSERTION_LEVEL_NORMAL)
                if (derived_state == InvalidAllelicState) {
                    derived_state = this_mutations_state;
                } else {
                    std::cout << "derived_state: " << derived_state << std::endl;
                    std::cout << "this_mutations_state: " << this_mutations_state << std::endl;
                    std::cout << "ancestral_state: " << ancestral_state << std::endl;
                    KASSERT(
                        (derived_state == this_mutations_state || derived_state == ancestral_state),
                        "Site is multiallelic. Allele frequency iterator only supports biallelic sites.",
                        tdt::assert::normal
                    );
                }
#endif

                AllelicState parent_state;
                // TODO this should be encapsulated in the GenomicSequenceStorage
                if (mutation.parent_mutation_id() == TSK_NULL) {
                    parent_state = ancestral_state;
                } else {
                    parent_state =
                        _sequence.mutation_by_id(asserting_cast<size_t>(mutation.parent_mutation_id())).allelic_state();
                }

                // TODO Make this branchless
                if (this_mutations_state != parent_state) {
                    if (this_mutations_state == ancestral_state) {
                        KASSERT(
                            num_samples_in_ancestral_state + num_samples_below_this_mutation <= num_samples(),
                            "There should never be more samples in the ancestral state than total samples.",
                            tdt::assert::light
                        );
                        num_samples_in_ancestral_state += num_samples_below_this_mutation;
                    } else {
                        KASSERT(
                            num_samples_in_ancestral_state >= num_samples_below_this_mutation,
                            "There should never be more samples in the derived state than total samples.",
                            tdt::assert::light
                        );
                        num_samples_in_ancestral_state -= num_samples_below_this_mutation;
                    }
                }
            }

            _frequency = num_samples_in_ancestral_state;
        }
    };

private:
    CompressedForest&             _forest;
    GenomicSequenceStorage const& _sequence;
};

// TODO encapsulate CompresedForest and GenomicSequence functions in this class
class SequenceForest {
public:
    // TODO Rename to GenomicSequenceStore
    SequenceForest(
        TSKitTreeSequence&&      tree_sequence,
        CompressedForest&&       compressed_forest,
        GenomicSequenceStorage&& genomic_sequence
    )
        : _tree_sequence(std::move(tree_sequence)),
          _forest(compressed_forest),
          _sequence(genomic_sequence) {}

    SequenceForest(tsk_treeseq_t&& ts_tree_sequence) : SequenceForest(TSKitTreeSequence(ts_tree_sequence)) {}

    SequenceForest(TSKitTreeSequence&& tree_sequence) : _tree_sequence(std::move(tree_sequence)) {
        ForestCompressor              forest_compressor(_tree_sequence);
        GenomicSequenceStorageFactory sequence_factory(_tree_sequence);
        _forest   = forest_compressor.compress(sequence_factory);
        _sequence = sequence_factory.move_storage();
    }

    AlleleFrequencies allele_frequencies() {
        return AlleleFrequencies(_forest, _sequence);
    }

    // TODO Make this const
    [[nodiscard]] double diversity() {
        auto const n  = num_samples();
        double     pi = 0.0;
        for (uint64_t frequency: allele_frequencies()) {
            pi += static_cast<double>(frequency * (n - frequency));
        }

        return pi / (0.5 * static_cast<double>(n * (n - 1)));
    }

    // TODO Make this const
    [[nodiscard]] uint64_t num_segregating_sites() {
        size_t num_segregating_sites = 0;
        for (uint64_t frequency: allele_frequencies()) {
            if (frequency < num_samples()) {
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
        // TODO are there formulas for computing these for efficiently?
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

    [[nodiscard]] GenomicSequenceStorage const& sequence() const {
        return _sequence;
    }

    [[nodiscard]] CompressedForest const& forest() const {
        return _forest;
    }

    [[nodiscard]] TSKitTreeSequence const& tree_sequence() const {
        return _tree_sequence;
    }

private:
    TSKitTreeSequence      _tree_sequence;
    CompressedForest       _forest;
    GenomicSequenceStorage _sequence;
};
