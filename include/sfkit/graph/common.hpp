#pragma once

#include <cstddef>

#include <kassert/kassert.hpp>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/utils/xxhash.hpp"

using NodeId = uint32_t;
using TreeId = uint32_t;
using EdgeId = uint32_t;

constexpr auto NodeId_bitwidth = sizeof(NodeId) * 8;
constexpr auto TreeId_bitwidth = sizeof(TreeId) * 8;
constexpr auto EdgeId_bitwidth = sizeof(EdgeId) * 8;

constexpr NodeId INVALID_NODE_ID = static_cast<NodeId>(-1);

class Edge {
public:
    Edge() = default;

    Edge(NodeId from, NodeId to) : _from(from), _to(to) {}

    NodeId from() const {
        KASSERT(_from != INVALID_NODE_ID, "Invalid node ID.", sfkit::assert::light);
        return _from;
    }

    void from(NodeId const from) {
        KASSERT(_from != INVALID_NODE_ID, "Invalid node ID.", sfkit::assert::light);
        _from = from;
    }

    NodeId to() const {
        KASSERT(_to != INVALID_NODE_ID, "Invalid node ID.", sfkit::assert::light);
        return _to;
    }

    void to(NodeId const to) {
        KASSERT(_to != INVALID_NODE_ID, "Invalid node ID.", sfkit::assert::light);
        _to = to;
    }

    template <class Archive>
    void serialize(Archive& ar) {
        ar(_from, _to);
    }

    bool operator==(Edge const& other) const noexcept {
        return _from == other._from && _to == other._to;
    }

private:
    NodeId _from = INVALID_NODE_ID;
    NodeId _to   = INVALID_NODE_ID;
};

struct EdgeHash {
    size_t operator()(Edge const& edge) const {
        return xxhash64(edge.to(), xxhash64(edge.from()));
    }
};

enum class TraversalOrder { Preorder, Postorder, Levelorder, Inorder, Unordered };
