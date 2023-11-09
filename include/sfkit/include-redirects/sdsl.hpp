#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wshadow"
#if defined(__GNUC__) && !defined(__clang__)
    #pragma GCC diagnostic ignored "-Wsuggest-override"
#endif
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#pragma GCC diagnostic pop
