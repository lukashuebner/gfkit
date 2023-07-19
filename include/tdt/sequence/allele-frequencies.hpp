#pragma once

#include <vector>

#include "tdt/graph/compressed-forest.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"

// TODO This class shares functionality with the AlleleFrequencySpectrum class. We should refactor this.
// The historical reason for this duplication is that AlleleFrequencySpectrum was written under the assumption of
// multiallelicity and this class was written under the assumption of biallelicity when I implemented the diversity
// statistics.
class AlleleFrequencies {
public:
    AlleleFrequencies(
        CompressedForest& compressed_forest, GenomicSequence const& sequence_store, SampleSet const& sample_set
    )
        : _forest(compressed_forest),
          _sequence(sequence_store),
          _num_samples_below(compressed_forest.compute_num_samples_below(sample_set)) {}

    auto begin() const {
        return allele_frequency_iterator{*this};
    }

    auto end() const {
        return allele_frequency_iterator::sentinel{};
    }

    auto cbegin() const {
        return allele_frequency_iterator{*this};
    }

    auto cend() const {
        return allele_frequency_iterator::sentinel{};
    }

    CompressedForest& forest() {
        return _forest;
    }

    GenomicSequence const& sequence() const {
        return _sequence;
    }
    
    NumSamplesBelow const& num_samples_below() const {
        return _num_samples_below;
    }

    class allele_frequency_iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = uint64_t;
        using pointer           = value_type*;
        using reference         = value_type&;
        struct sentinel {};

        allele_frequency_iterator(AlleleFrequencies const& freqs) : _freqs(freqs), _site(0) {
            _update_state();
        }

        // The num_samples below should not be recomputed every time we copy the iterator

        [[nodiscard]] SampleId num_samples() const {
            return _freqs._num_samples_below.num_samples_in_sample_set();
        }

        [[nodiscard]] SiteId num_sites() const {
            return _freqs._sequence.num_sites();
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

        [[nodiscard]] bool operator==(allele_frequency_iterator const& other) const {
            return &_freqs == &other._freqs && _site == other._site;
        }

        [[nodiscard]] bool operator!=(allele_frequency_iterator const& other) const {
            return !(*this == other);
        }

        [[nodiscard]] bool operator==(sentinel) const {
            return _site == num_sites();
        }

        reference operator*() {
            return _frequency;
        }

        pointer operator->() {
            return &_frequency;
        }

    private:
        AlleleFrequencies const& _freqs;
        uint64_t                 _frequency;
        SiteId                   _site;

        void _update_state() {
            // TODO Generalize this to multiallelic sites
            KASSERT(_site < num_sites(), "Site index out of bounds", tdt::assert::light);

            uint64_t           num_samples_in_ancestral_state = num_samples();
            AllelicState const ancestral_state                = _freqs._sequence.ancestral_state(_site);

#if KASSERT_ENABLED(TDT_ASSERTION_LEVEL_NORMAL)
            AllelicState derived_state = InvalidAllelicState;
#endif

            for (auto const& mutation: _freqs._sequence.mutations_at_site(_site)) {
                auto const num_samples_below_this_mutation = _freqs._num_samples_below(mutation.subtree_id());
                auto const this_mutations_state            = mutation.allelic_state();

#if KASSERT_ENABLED(TDT_ASSERTION_LEVEL_NORMAL)
                if (derived_state == InvalidAllelicState) {
                    derived_state = this_mutations_state;
                } else {
                    // std::cout << "derived_state: " << derived_state << std::endl;
                    // std::cout << "this_mutations_state: " << this_mutations_state << std::endl;
                    // std::cout << "ancestral_state: " << ancestral_state << std::endl;
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
                        _freqs._sequence.mutation_by_id(asserting_cast<size_t>(mutation.parent_mutation_id()))
                            .allelic_state();
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
    CompressedForest&      _forest;
    GenomicSequence const& _sequence;
    NumSamplesBelow        _num_samples_below;
};
