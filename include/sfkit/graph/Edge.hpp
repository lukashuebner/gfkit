#pragma once

#include <kassert/kassert.hpp>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/graph/primitives.hpp"
#include "sfkit/utils/xxhash.hpp"

namespace sfkit::graph {
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
        return sfkit::utils::xxhash64(edge.to(), sfkit::utils::xxhash64(edge.from()));
    }
};
} // namespace sfkit::graph
