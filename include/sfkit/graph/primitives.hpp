#pragma once

#include <cstddef>
#include <cstdint>

namespace sfkit::graph {
using NodeId = uint32_t;
using TreeId = uint32_t;
using EdgeId = uint32_t;

constexpr auto NodeId_bitwidth = sizeof(NodeId) * 8;
constexpr auto TreeId_bitwidth = sizeof(TreeId) * 8;
constexpr auto EdgeId_bitwidth = sizeof(EdgeId) * 8;

constexpr NodeId INVALID_NODE_ID = static_cast<NodeId>(-1);

enum class TraversalOrder { Preorder, Postorder, Levelorder, Inorder, Unordered };
} // namespace sfkit::graph
