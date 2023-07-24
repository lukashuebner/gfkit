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

class PerfectNumericHasher {
public:
    using Idx = uint8_t;

    static inline Idx to_idx(AllelicState state) {
#if KASSERT_ENABLED(TDT_ASSERTION_LEVEL_LIGHT)
        KASSERT(
            state == '0' || state == '1' || state == '2' || state == '3' || state == '4' || state == '5' || state == '6'
                || state == '7' || state == '8' || state == '9',
            "Invalid numeric state: " << state,
            tdt::assert::light
        );
#endif

        char const idx = state - '0';
        KASSERT((idx >= 0 && idx < num_states), "Invalid state: " << state, tdt::assert::light);
        return static_cast<Idx>(idx);
    }

    static constexpr Idx num_states = 4;
};

class PerfectDNAHasher {
public:
    using Idx = uint8_t;

    static inline Idx to_idx(AllelicState state) {
#if KASSERT_ENABLED(TDT_ASSERTION_LEVEL_LIGHT)
        KASSERT(
            state == 'A' || state == 'C' || state == 'G' || state == 'T',
            "Invalid DNA state: " << state,
            tdt::assert::light
        );
#endif
        // The characters A, C, T and G differ pairwise in the second and third least significant bit.
        // A 0100 0001
        // C 0100 0011
        // G 0100 0111
        // T 0101 0100
        //         ^^
        constexpr char ISOLATE_BITMASK = 0b0000'0110;
        constexpr Idx  SHIFT           = 1;
        Idx            idx             = asserting_cast<Idx>(state & ISOLATE_BITMASK) >> SHIFT;
        KASSERT(idx < num_states, "Invalid state: " << state, tdt::assert::light);
        return idx;
    }

    static constexpr Idx num_states = 4;
};
