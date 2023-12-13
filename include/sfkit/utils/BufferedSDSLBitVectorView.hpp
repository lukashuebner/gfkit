#pragma once

#include <kassert/kassert.hpp>
#include <sdsl/bit_vectors.hpp>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/utils/bitpacking.hpp"
#include "sfkit/utils/checking_casts.hpp"

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
            --_remaining_bits_in_buffer;

            return *this;
        }

        iterator operator++(int) {
            iterator tmp(*this);

            this->operator++();
            return tmp;
        }

        [[nodiscard]] bool operator==(sentinel const) const {
            return _is_end;
        }

        [[nodiscard]] bool operator!=(sentinel const) const {
            return !_is_end;
        }

        [[nodiscard]] reference operator*() {
            KASSERT(_remaining_bits_in_buffer < BUFFER_SIZE);
            return _buffer[_remaining_bits_in_buffer];
        }

        [[nodiscard]] pointer operator->() {
            KASSERT(_remaining_bits_in_buffer < BUFFER_SIZE);
            return &(_buffer[_remaining_bits_in_buffer]);
        }

    private:
        using buffer_size_t                        = uint8_t;
        static constexpr buffer_size_t BUFFER_SIZE = 8;
        using buffer_t                             = std::array<bool, BUFFER_SIZE>;

        // static constexpr buffer_size_t BUFFER_SIZE   = sizeof(buffer_t) * 8;
        // static constexpr buffer_t      RIGHTMOST_BIT = 1ULL;
        // static constexpr buffer_t      LEFTMOST_BIT  = 1ULL << (BUFFER_SIZE - 1);

        sdsl::bit_vector const&     _bit_vector;
        sdsl::bit_vector::size_type _bit_vector_idx;
        sdsl::bit_vector::size_type _bit_vector_size;
        buffer_t                    _buffer;
        buffer_size_t               _remaining_bits_in_buffer = 0;
        bool                        _is_end                   = false;

        void _refill_buffer() {
            KASSERT(_bit_vector_idx <= _bit_vector_size, "Bit vector index is out of bounds", sfkit::assert::light);
            buffer_size_t const bits_read =
                asserting_cast<buffer_size_t>(std::min<size_t>(BUFFER_SIZE, _bit_vector_size - _bit_vector_idx));
            if (bits_read == 0) [[unlikely]] {
                _is_end = true;
            } else {
                uint8_t bits = asserting_cast<uint8_t>(_bit_vector.get_int(_bit_vector_idx, bits_read));
                KASSERT(BUFFER_SIZE >= bits_read);
                bits = asserting_cast<buffer_size_t>(
                    bits << asserting_cast<buffer_size_t>(
                        static_cast<buffer_size_t>(BUFFER_SIZE) - static_cast<buffer_size_t>(bits_read)
                    )
                );
                unpack8bools(bits, _buffer.data());
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
