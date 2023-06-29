#pragma once

#include <kassert/kassert.hpp>
#include <tskit.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/checking_casts.hpp"

// TODO Use sfkit namespace

using SiteId = uint32_t;
// TODO Use only 2 bits per site
using AllelicState                         = char;
constexpr AllelicState InvalidAllelicState = -1;

// TODO Implement something safer here (e.g. what happens if '1' and 'A' gets mixed in the same dataset?)
class PerfectNumericHasher {
public:
    static inline unsigned char to_idx(AllelicState state) {
        char const idx = state - '0';
        KASSERT((idx >= 0 && idx < num_states), "Invalid state: " << state, tdt::assert::light);
        return static_cast<unsigned char>(idx);
    }

    static constexpr unsigned char num_states = 4;
};

class PerfectDNAHasher {
public:
    static inline unsigned char to_idx(AllelicState state) {
        // The characters A, C, T and G differ pairwise in the second and third least significant bit.
        // A 0100 0001
        // C 0100 0011
        // G 0100 0111
        // T 0101 0100
        //         ^^
        constexpr char          ISOLATE_BITMASK = 0b0000'0110;
        constexpr unsigned char SHIFT           = 1;
        unsigned char           idx             = asserting_cast<unsigned char>(state & ISOLATE_BITMASK) >> SHIFT;
        KASSERT(idx < num_states, "Invalid state: " << state, tdt::assert::light);
        return idx;
    }

    static constexpr unsigned int num_states = 4;
};
