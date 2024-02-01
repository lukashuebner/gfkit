#pragma once

#include <array>
#include <cstdint>
#include <numeric>
#include <vector>

#include <kassert/kassert.hpp>
#include <tskit.h>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/graph/primitives.hpp"
#include "sfkit/samples/primitives.hpp"
#include "sfkit/utils/checking_casts.hpp"
#include "sfkit/utils/concepts.hpp"

namespace sfkit::samples {

using sfkit::utils::asserting_cast;
using sfkit::utils::IterableInput;

using SampleSetId = uint8_t;

class SampleSet {
public:
    class const_iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = SampleId;
        using pointer           = value_type*;
        using reference         = value_type&;
        struct sentinel {};

        using underlying_iterator = std::vector<bool>::const_iterator;

        const_iterator(underlying_iterator begin, underlying_iterator end) : _it(begin), _end(end), _sample_id(0) {
            _forward_to_next_contained_sample();
        }

        const_iterator& operator++() {
            ++_it;
            ++_sample_id;
            _forward_to_next_contained_sample();
            return *this;
        }

        const_iterator operator++(int) {
            const_iterator tmp(*this);
            this->         operator++();
            return tmp;
        }

        bool operator==(sentinel) {
            return _it == _end;
        }

        reference operator*() {
            return _sample_id;
        }

        pointer operator->() {
            return &_sample_id;
        }

    private:
        underlying_iterator _it;
        underlying_iterator _end;
        SampleId            _sample_id;

        void _forward_to_next_contained_sample() {
            while (_it != _end && !*_it) {
                _it++;
                _sample_id++;
            }
        }
    };

    SampleSet(SampleId const num_samples_in_dag) {
        _samples.resize(num_samples_in_dag);
    };

    const_iterator begin() const {
        return const_iterator(_samples.cbegin(), _samples.cend());
    }

    const_iterator::sentinel end() const {
        return const_iterator::sentinel();
    }

    SampleSet& add(SampleId const sample_id) {
        KASSERT(sample_id < _samples.size(), "Sample ID out of bounds.", sfkit::assert::light);
        _samples[sample_id] = true;
        return *this;
    }

    template <IterableInput IterableInput>
    SampleSet& add(IterableInput const& input) {
        for (auto&& sample_id: input) {
            this->add(sample_id);
        }
        return *this;
    }

    SampleSet& remove(SampleId const sample_id) {
        KASSERT(sample_id < _samples.size(), "Sample ID out of bounds.", sfkit::assert::light);
        _samples[sample_id] = false;
        return *this;
    }

    template <IterableInput IterableInput>
    SampleSet& remove(IterableInput const& input) {
        for (auto&& sample_id: input) {
            this->remove(sample_id);
        }
        return *this;
    }

    SampleSet& clear() {
        std::fill(_samples.begin(), _samples.end(), false);
        return *this;
    }

    [[nodiscard]] SampleId operator[](SampleId const sample_id) const {
        KASSERT(sample_id < _samples.size(), "Sample ID out of bounds.", sfkit::assert::light);
        return _samples[sample_id];
    }

    [[nodiscard]] SampleId overall_num_samples() const {
        return asserting_cast<SampleId>(_samples.size());
    }

    [[nodiscard]] SampleId popcount() const {
        auto const popcount = std::reduce(_samples.begin(), _samples.end(), 0);
        KASSERT(
            std::cmp_less_equal(popcount, _samples.size()),
            "The popcount is larger than the number of samples.",
            sfkit::assert::light
        );
        KASSERT(std::cmp_greater_equal(popcount, 0), "The popcount is negative.", sfkit::assert::light);
        return asserting_cast<SampleId>(popcount);
    }

    [[nodiscard]] SampleSet inverse() const {
        SampleSet inverse(overall_num_samples());
        for (SampleId sample_id = 0; sample_id < _samples.size(); sample_id++) {
            if (!(*this)[sample_id]) {
                inverse.add(sample_id);
            }
        }
        return inverse;
    }

    [[nodiscard]] std::vector<tsk_id_t> to_tsk_samples() const {
        std::vector<tsk_id_t> tsk_samples;

        tsk_samples.reserve(popcount());
        for (auto&& sample_id: *this) {
            tsk_samples.emplace_back(asserting_cast<tsk_id_t>(sample_id));
        }

        KASSERT(
            tsk_samples.size() == popcount(),
            "The number of samples in the TSK samples vector does not match the popcount.",
            sfkit::assert::light
        );

        return tsk_samples;
    }

private:
    std::vector<bool> _samples;
};

template <size_t N>
using SetOfSampleSets = std::array<std::reference_wrapper<SampleSet const>, N>;

} // namespace sfkit::samples
