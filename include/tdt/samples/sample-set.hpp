#pragma once

#include <cstdint>
#include <numeric>
#include <vector>

#include <kassert/kassert.hpp>
#include <tskit.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/checking_casts.hpp"
#include "tdt/utils/concepts.hpp"

using SampleId = uint32_t;

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
            _it++;
            _sample_id++;
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

        bool operator==(const_iterator const& other) const {
            return _it == other._it && _end == other._end;
        }

        bool operator!=(const_iterator const& other) const {
            return !(*this == other);
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

    void add(SampleId const sample_id) {
        KASSERT(sample_id < _samples.size(), "Sample ID out of bounds.", tdt::assert::light);
        _samples[sample_id] = true;
    }

    template <IterableInput IterableInput>
    void add(IterableInput input) {
        for (auto&& sample_id: input) {
            this->add(sample_id);
        }
    }

    void remove(SampleId const sample_id) {
        KASSERT(sample_id < _samples.size(), "Sample ID out of bounds.", tdt::assert::light);
        _samples[sample_id] = false;
    }

    template <IterableInput IterableInput>
    void remove(IterableInput input) {
        for (auto&& sample_id: input) {
            this->remove(sample_id);
        }
    }

    void clear() {
        std::fill(_samples.begin(), _samples.end(), false);
    }

    SampleId operator[](SampleId const sample_id) const {
        KASSERT(sample_id < _samples.size(), "Sample ID out of bounds.", tdt::assert::light);
        return _samples[sample_id];
    }

    SampleId num_nodes_in_dag() const {
        return asserting_cast<SampleId>(_samples.size());
    }

    SampleId popcount() const {
        return asserting_cast<SampleId>(std::accumulate(_samples.begin(), _samples.end(), 0));
    }

    // Building the inverse is not that trivial, as we don't know which nodes are inner nodes and which are leaves.
    // SampleSet build_inverse() const {
    //     // Let's see how well the compiler optimizes this and optimize if the profiler tells us that this is too slow.
    //     SampleSet inverse(overall_num_samples());
    //     for (SampleId sample_id = 0; sample_id < _samples.size(); ++sample_id) {
    //         if (!_samples[sample_id]) {
    //             inverse.add(sample_id);
    //         }
    //     }

    //     return inverse;
    // }

    std::vector<tsk_id_t> as_tsk_id_t_vector() const {
        std::vector<tsk_id_t> result;
        result.reserve(this->popcount());
        for (SampleId sample_id = 0; sample_id < _samples.size(); ++sample_id) {
            if (_samples[sample_id]) {
                result.push_back(asserting_cast<tsk_id_t>(sample_id));
            }
        }
        return result;
    }

private:
    // TODO Profile and check if a compressed bitset would be faster?
    std::vector<bool> _samples;
};
