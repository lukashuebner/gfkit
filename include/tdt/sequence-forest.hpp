
#include <kassert/kassert.hpp>
#include <tskit.h>

#include "assertion_levels.hpp"
#include "graph/compressed-forest.hpp"
#include "sequence/genomic-sequence-storage.hpp"

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
        size_t                        _site;
        uint64_t                      _frequency;

        void _update_state() {
            // TODO Generalize this to multiallelic sites
            KASSERT(_site < num_sites(), "Site index out of bounds", tdt::assert::light);

            uint64_t           num_samples_in_ancestral_state = num_samples();
            AllelicState const ancestral_state                = _sequence.ancestral_state(_site);

            // TODO KASSERT() for biallelicity
            for (auto const& mutation: _sequence.mutations_at_site(_site)) {
                auto const num_samples_below_this_mutation = _forest.num_samples_below(mutation.subtree_id());
                auto const this_mutations_state            = mutation.allelic_state();

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
                            num_samples_in_ancestral_state + num_samples_below_this_mutation < num_samples(),
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
    SequenceForest(CompressedForest& compressed_forest, GenomicSequenceStorage& genomic_sequence)
        : _forest(compressed_forest),
          _sequence(genomic_sequence) {}

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

    [[nodiscard]] uint64_t num_sites() const {
        return _sequence.num_sites();
    }

    [[nodiscard]] uint64_t num_samples() const {
        return _forest.num_samples();
    }

private:
    CompressedForest&       _forest;
    GenomicSequenceStorage& _sequence;
};
