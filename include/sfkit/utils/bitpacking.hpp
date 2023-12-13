#pragma once

#include <cstdint>
#include <cstring>

// https://stackoverflow.com/questions/8461126/how-to-create-a-byte-out-of-8-bool-values-and-vice-versa/51750902#51750902
inline void unpack8bools(uint8_t bits, bool* dest) {
    // on little-endian,  a[0] = (b>>7) & 1  like printing order
    auto     MAGIC = 0x8040201008040201ULL; // for opposite order, byte-reverse this
    auto     MASK  = 0x8080808080808080ULL;
    uint64_t t     = ((MAGIC * bits) & MASK) >> 7;
    memcpy(dest, &t, sizeof t); // store 8 bytes without UB
}
