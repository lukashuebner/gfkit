#pragma once

#include <cstddef>

using NodeId = std::size_t;
using SubtreeId = NodeId;
using TreeId = NodeId;

class Edge {
public:
    Edge(NodeId from, NodeId to) : _from(from), _to(to) {}

    NodeId from() const {
        return _from;
    }

    void from(NodeId const from) {
        _from = from;
    }

    NodeId to() const {
        return _to;
    }

    void to(NodeId const to) {
        _to = to;
    }

private:
    NodeId _from;
    NodeId _to;
};

enum class TraversalOrder { Preorder, Postorder, Levelorder, Inorder, Unordered };
