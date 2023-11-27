#pragma once

#include <variant>
#include <vector>

// TODO Generalize to be able to use DAG and BP based forests
#include "sfkit/dag/DAGCompressedForest.hpp"
#include "sfkit/samples/NumSamplesBelowAccessor.hpp"
#include "sfkit/samples/NumSamplesBelowFactory.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"
#include "sfkit/utils/always_false_v.hpp"

namespace sfkit::sequence {

using namespace sfkit::samples;
using sfkit::dag::DAGCompressedForest;
using sfkit::samples::DAGNumSamplesBelow;
using sfkit::samples::NumSamplesBelowAccessor;
using sfkit::utils::always_false_v;

template <
    typename AllelicStatePerfectHasher                = PerfectDNAHasher,
    NumSamplesBelowAccessorC NumSamplesBelowAccessorT = NumSamplesBelowAccessor<DAGNumSamplesBelow<1>>>
class AlleleFrequencies {
public:
    class BiallelicFrequency {
    public:
        BiallelicFrequency(SampleId num_ancestral) noexcept : _num_ancestral(num_ancestral) {}

        [[nodiscard]] SampleId num_ancestral() const {
            return _num_ancestral;
        }

        // Always compare the number of ancestral samples. Only compare the ancestral state and derived state if both
        // values have it set.
        bool operator==(BiallelicFrequency const& other) const {
            return _num_ancestral == other._num_ancestral;
        }

    private:
        SampleId _num_ancestral; // There can never be more derived samples than total samples.
    };

    class MultiallelicFrequency {
    public:
        using iterator       = typename std::array<SampleId, AllelicStatePerfectHasher::num_states>::iterator;
        using const_iterator = typename std::array<SampleId, AllelicStatePerfectHasher::num_states>::const_iterator;

        using Idx = typename AllelicStatePerfectHasher::Idx;

        MultiallelicFrequency(AllelicState ancestral_state = InvalidAllelicState) : _ancestral_state(ancestral_state) {
            _num_samples_in_state.fill(0);
        }

        template <typename... T>
        MultiallelicFrequency(AllelicState ancestral_state, T... ts)
            : _num_samples_in_state{ts...},
              _ancestral_state(ancestral_state) {}

        [[nodiscard]] SampleId operator[](Idx idx) const {
            KASSERT(idx < AllelicStatePerfectHasher::num_states, "Index out of bounds", sfkit::assert::light);
            return _num_samples_in_state[idx];
        }

        SampleId& operator[](Idx idx) {
            KASSERT(idx < AllelicStatePerfectHasher::num_states, "Index out of bounds", sfkit::assert::light);
            return _num_samples_in_state[idx];
        }

        iterator begin() {
            return _num_samples_in_state.begin();
        }

        iterator end() {
            return _num_samples_in_state.end();
        }

        const_iterator begin() const noexcept {
            return _num_samples_in_state.begin();
        }

        const_iterator end() const noexcept {
            return _num_samples_in_state.end();
        }

        [[nodiscard]] bool valid(SampleId expected_num_samples) const {
            auto const sum = std::accumulate(
                _num_samples_in_state.begin(),
                _num_samples_in_state.end(),
                SampleId(0),
                [](SampleId running_sum, SampleId num_samples_in_state) { return running_sum + num_samples_in_state; }
            );
            return sum == expected_num_samples;
        }

        [[nodiscard]] AllelicState ancestral_state() const {
            KASSERT(
                _ancestral_state != InvalidAllelicState,
                "The ancestral state has not been set.",
                sfkit::assert::light
            );
            return _ancestral_state;
        }

        [[nodiscard]] Idx ancestral_state_idx() const {
            KASSERT(
                _ancestral_state != InvalidAllelicState,
                "The ancestral state has not been set.",
                sfkit::assert::light
            );
            return AllelicStatePerfectHasher::to_idx(_ancestral_state);
        }

        // Always compare the number of samples per state. Only compare the ancestral state if both values have it set.
        bool operator==(MultiallelicFrequency const& other) const {
            return _num_samples_in_state == other._num_samples_in_state
                   && (_ancestral_state == InvalidAllelicState || other._ancestral_state == InvalidAllelicState
                       || _ancestral_state == other._ancestral_state);
        }

        static constexpr auto num_states = AllelicStatePerfectHasher::num_states;

    private:
        using AllelicStateFrequencies = std::array<SampleId, AllelicStatePerfectHasher::num_states>;
        AllelicStateFrequencies _num_samples_in_state;
        AllelicState            _ancestral_state;
    };

    using AlleleFrequency = std::variant<BiallelicFrequency, MultiallelicFrequency>;

    AlleleFrequencies(
        DAGCompressedForest& compressed_forest, GenomicSequence const& sequence_store, SampleSet const& sample_set
    )
        : _forest(compressed_forest),
          _sequence(sequence_store),
          _num_samples_below(NumSamplesBelowFactory::build(compressed_forest.postorder_edges(), sample_set)) {}

    AlleleFrequencies(
        DAGCompressedForest&            compressed_forest,
        GenomicSequence const&          sequence_store,
        NumSamplesBelowAccessorT const& num_samples_below
    )
        : _forest(compressed_forest),
          _sequence(sequence_store),
          _num_samples_below(num_samples_below) {}

    [[nodiscard]] auto begin() const {
        return allele_frequency_iterator{*this};
    }

    [[nodiscard]] auto end() const {
        return typename allele_frequency_iterator::sentinel{};
    }

    [[nodiscard]] auto cbegin() const {
        return allele_frequency_iterator{*this};
    }

    [[nodiscard]] auto cend() const {
        return typename allele_frequency_iterator::sentinel{};
    }

    [[nodiscard]] DAGCompressedForest& forest() {
        return _forest;
    }

    [[nodiscard]] GenomicSequence const& sequence() const {
        return _sequence;
    }

    [[nodiscard]] NumSamplesBelowAccessorT const& num_samples_below() const {
        return _num_samples_below;
    }

    [[nodiscard]] SampleId num_samples_in_sample_set() const {
        return _num_samples_below.num_samples_in_sample_set();
    }

    template <typename BiallelicVisitor, typename MultiallelicVisitor>
    void visit(BiallelicVisitor biallelic_vistor, MultiallelicVisitor multiallelic_visitor) const noexcept {
        for (AlleleFrequency const& this_sites_state: *this) {
            std::visit(
                // TODO This noexcept should be conditional on the noexcept of the visitors.
                [biallelic_vistor, multiallelic_visitor](auto const& arg) noexcept -> void {
                    using T = std::decay_t<decltype(arg)>;
                    if constexpr (std::is_same_v<T, BiallelicFrequency>) {
                        biallelic_vistor(arg);
                    } else if constexpr (std::is_same_v<T, MultiallelicFrequency>) {
                        multiallelic_visitor(arg);
                    } else {
                        static_assert(
                            always_false_v<T>,
                            "Returned type of allele_frequencies is neither bi- nor multiallelic."
                        );
                    }
                },
                this_sites_state
            );
        }
    }

    class allele_frequency_iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = AlleleFrequency;
        using pointer           = value_type*;
        using reference         = value_type&;
        struct sentinel {};

        allele_frequency_iterator(AlleleFrequencies const& freqs)
            : _freqs(freqs),
              _state(BiallelicFrequency(0)),
              _site(0) {
            _update_state();
        }

        [[nodiscard]] SampleId num_samples_in_sample_set() const {
            return _freqs._num_samples_below.num_samples_in_sample_set();
        }

        [[nodiscard]] SiteId num_sites() const {
            return _freqs._sequence.num_sites();
        }

        // TODO Enable iteration only over sites with mutations
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
            return _state;
        }

        pointer operator->() {
            return &_state;
        }

        void force_multiallelicity() {
            if (std::holds_alternative<BiallelicFrequency>(_state)) {
                AllelicState const ancestral_state   = _freqs._sequence.ancestral_state(_site);
                auto const         mutations_at_site = _freqs._sequence.mutations_at_site(_site);
                _update_state_multiallelic(ancestral_state, mutations_at_site);
            }
        }

    private:
        AlleleFrequencies const& _freqs;
        AlleleFrequency          _state;
        SiteId                   _site;

        // Maybe it's faster (it's certainly simpler) to just use multiallelic iterator for all sites.
        void _update_state() {
            // TODO Can we improve performance by knowing if a site is bi- or multiallelic?
            KASSERT(_site < num_sites(), "Site index out of bounds", sfkit::assert::light);

            SampleId           num_ancestral     = num_samples_in_sample_set();
            AllelicState const ancestral_state   = _freqs._sequence.ancestral_state(_site);
            auto const         mutations_at_site = _freqs._sequence.mutations_at_site(_site);

            // Handle case where there are no mutations at this site
            // TODO Remove this once we're iterating only over sites with mutations
            if (mutations_at_site.empty()) {
                _state = BiallelicFrequency(num_ancestral);
                return;
            } else {
                auto         mutation_it   = mutations_at_site.begin();
                AllelicState derived_state = InvalidAllelicState;
                // TODO Filter those mutations during construction of the compressed forest
                do { // The first mutations could be from the ancestral state to the ancestral state
                    derived_state = mutation_it->allelic_state();
                } while (derived_state == ancestral_state && ++mutation_it != mutations_at_site.end());

                while (mutation_it != mutations_at_site.end()) {
                    auto const mutation                        = *mutation_it;
                    auto const num_samples_below_this_mutation = _freqs._num_samples_below(mutation.node_id());
                    auto const this_mutations_state            = mutation.allelic_state();

                    // If this mutation is towards neither the derived nor the ancestral state but there is at
                    // least one sample in this sample set below it, we have to compute the state using the
                    // algorithm supporting multiallelicity.
                    if (this_mutations_state != derived_state && this_mutations_state != ancestral_state
                        && num_samples_below_this_mutation > 0) [[unlikely]] {
                        _update_state_multiallelic(ancestral_state, mutations_at_site);
                        return;
                    }

                    // TODO Check if making this branchless is faster.
                    if (this_mutations_state != mutation.parent_state()) {
                        if (this_mutations_state == ancestral_state) {
                            KASSERT(
                                num_ancestral + num_samples_below_this_mutation <= num_samples_in_sample_set(),
                                "There should never be more samples in the ancestral state than total samples.",
                                sfkit::assert::light
                            );
                            num_ancestral += num_samples_below_this_mutation;
                        } else {
                            KASSERT(
                                num_ancestral >= num_samples_below_this_mutation,
                                "There should never be more samples in the derived state than total samples.",
                                sfkit::assert::light
                            );
                            num_ancestral -= num_samples_below_this_mutation;
                        }
                    }
                    mutation_it++;
                }

                _state = BiallelicFrequency(num_ancestral);
                return;
            }
        }

        void _update_state_multiallelic(AllelicState ancestral_state, MutationView const& mutations_at_site) {
            // For this site, compute how many of the samples have which derived or ancestral state.
            MultiallelicFrequency state_freqs(ancestral_state);

            // Before looking at any mutations, all samples are in the ancestral state
            auto const idx_of_ancestral_state   = AllelicStatePerfectHasher::to_idx(ancestral_state);
            state_freqs[idx_of_ancestral_state] = num_samples_in_sample_set();

            for (auto const& mutation: mutations_at_site) {
                auto const num_samples_below_this_mutation = _freqs._num_samples_below(mutation.node_id());
                auto const idx_of_mutations_state = AllelicStatePerfectHasher::to_idx(mutation.allelic_state());
                auto const idx_of_parents_state   = AllelicStatePerfectHasher::to_idx(mutation.parent_state());

                state_freqs[idx_of_mutations_state] += num_samples_below_this_mutation;
                state_freqs[idx_of_parents_state] -= num_samples_below_this_mutation;
                KASSERT(
                    state_freqs.valid(num_samples_in_sample_set()),
                    "The number of samples per state does not sum up to the total number of samples.",
                    sfkit::assert::normal
                );
            }

            _state = state_freqs;
        }
    };

private:
    DAGCompressedForest&     _forest;
    GenomicSequence const&   _sequence;
    NumSamplesBelowAccessorT _num_samples_below;
};
} // namespace sfkit::sequence
