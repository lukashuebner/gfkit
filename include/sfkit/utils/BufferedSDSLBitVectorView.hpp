#pragma once

#include <kassert/kassert.hpp>
#include <sdsl/bit_vectors.hpp>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/checking_casts.hpp"

// TODO Use namespaces for all of sfkit
namespace sfkit::utils {
class BufferedSDSLBitVectorView {
public:
    BufferedSDSLBitVectorView(sdsl::bit_vector const& bit_vector) : _bit_vector(bit_vector) {}

    class iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = bool;
        using pointer           = value_type*;
        using reference         = value_type&;
        struct sentinel {};

        iterator(sdsl::bit_vector const& bit_vector)
            : _bit_vector(bit_vector),
              _bit_vector_idx(0),
              _bit_vector_size(bit_vector.size()) {
            _refill_buffer();
            ++*this; // Extract the first bit
        }

        iterator& operator++() {
            if (_remaining_bits_in_buffer == 0) [[unlikely]] {
                _refill_buffer();
            }
            _value = _buffer & RIGHTMOST_BIT;
            _buffer >>= 1;
            --_remaining_bits_in_buffer;

            return *this;
        }

        iterator operator++(int) {
            iterator tmp(*this);

            this->operator++();
            return tmp;
        }

        // [[nodiscard]] bool operator==(iterator const& other) const {
        //     return _bit_vector_idx == other._bit_vector_idx
        //            && _remaining_bits_in_buffer == other._remaining_bits_in_buffer && _buffer == other._buffer;
        // }

        // [[nodiscard]] bool operator!=(iterator const& other) const {
        //     return !(*this == other);
        // }

        [[nodiscard]] bool operator==(sentinel const) {
            return _is_end;
        }

        [[nodiscard]] reference operator*() {
            return _value;
        }

        [[nodiscard]] pointer operator->() {
            return &_value;
        }

    private:
        using buffer_t      = uint64_t;
        using buffer_size_t = uint8_t;

        static constexpr buffer_size_t BUFFER_SIZE   = sizeof(buffer_t) * 8;
        static constexpr buffer_t      RIGHTMOST_BIT = 1ULL;
        static constexpr buffer_t      LEFTMOST_BIT  = 1ULL << (BUFFER_SIZE - 1);

        bool                        _value;
        sdsl::bit_vector const&     _bit_vector;
        sdsl::bit_vector::size_type _bit_vector_idx;
        sdsl::bit_vector::size_type _bit_vector_size;
        buffer_t                    _buffer;
        buffer_size_t               _remaining_bits_in_buffer = 0;
        bool                        _is_end                   = false;

        void _refill_buffer() {
            // TODO What happens if there aren't enough bits left in the bit_vector?
            KASSERT(_bit_vector_idx <= _bit_vector_size, "Bit vector index is out of bounds", sfkit::assert::light);
            auto const bits_read =
                asserting_cast<decltype(BUFFER_SIZE)>(std::min<size_t>(BUFFER_SIZE, _bit_vector_size - _bit_vector_idx)
                );
            if (bits_read == 0) [[unlikely]] {
                _is_end = true;
            } else {
                _buffer                   = _bit_vector.get_int(_bit_vector_idx, bits_read);
                _remaining_bits_in_buffer = bits_read;
                _bit_vector_idx += bits_read;
            }
        }
    };

    iterator begin() const {
        return iterator(_bit_vector);
    }

    iterator::sentinel end() const {
        return iterator::sentinel();
    }

private:
    sdsl::bit_vector const& _bit_vector;
};
} // namespace sfkit::utils
