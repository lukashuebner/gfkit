#pragma once

#include <span>

#include "sfkit/sequence/Mutation.hpp"

namespace sfkit::tskit {

// TODO Rename to MutationView to be in line with the naming of the other views
using TskMutationView = std::span<tsk_mutation_t const>;

} // namespace sfkit::tskit
