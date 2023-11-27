#pragma once

#include <iostream>
#include <vector>

#include "sfkit/utils/checking_casts.hpp"

// TODO rename this folder to io/ and introduce subfolder

namespace sfkit::io::utils {

using sfkit::utils::asserting_cast;

template <typename POD>
std::ostream& serialize(std::ostream& os, std::vector<POD> const& v) {
    // this only works on built in data types (PODs)
    static_assert(
        std::is_trivial_v<POD> && std::is_standard_layout_v<POD>,
        "Can only serialize POD types with this function"
    );

    auto size = v.size();
    os.write(reinterpret_cast<char const*>(&size), sizeof(size));
    os.write(reinterpret_cast<char const*>(v.data()), asserting_cast<std::streamsize>(v.size() * sizeof(POD)));
    return os;
}

template <typename POD>
std::istream& deserialize(std::istream& is, std::vector<POD>& v) {
    static_assert(
        std::is_trivial_v<POD> && std::is_standard_layout_v<POD>,
        "Can only deserialize POD types with this function"
    );

    decltype(v.size()) size;
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
    v.resize(size);
    is.read(reinterpret_cast<char*>(v.data()), asserting_cast<std::streamsize>(v.size() * sizeof(POD)));
    return is;
}
} // namespace sfkit::io::utils
