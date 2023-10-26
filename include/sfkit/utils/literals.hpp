#pragma once

inline unsigned char operator""_uc(unsigned long long value) {
    return static_cast<unsigned char>(value);
}
