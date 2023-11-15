#pragma once

#if defined(__GNUC__) && !defined(__clang__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wnoexcept"
#endif
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#if defined(__GNUC__) && !defined(__clang__)
    #pragma GCC diagnostic pop
#endif
