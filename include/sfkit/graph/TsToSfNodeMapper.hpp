#pragma once

#include <algorithm>
#include <optional>
#include <vector>

#include "sfkit/graph/primitives.hpp"
#include "sfkit/utils/concepts.hpp"

namespace sfkit::graph {
template <typename TsNodeToSubtreeMapper, typename SubtreeToSfNodeMapper>
class TsToSfNodeMapper {
public:
    TsToSfNodeMapper(
        TsNodeToSubtreeMapper const& ts_node_to_subtree_mapper, SubtreeToSfNodeMapper const& subtree_to_sf_node_mapper
    )
        : _ts_node_to_subtree_mapper(ts_node_to_subtree_mapper),
          _subtree_to_sf_node_mapper(subtree_to_sf_node_mapper) {}

    auto operator()(auto const& ts_node_id) const {
        return _subtree_to_sf_node_mapper[_ts_node_to_subtree_mapper[ts_node_id]];
    }

    std::optional<NodeId> try_map(auto const& ts_node_id) const {
        auto const subtree_id = _ts_node_to_subtree_mapper[ts_node_id];

        auto const sf_node_it = _subtree_to_sf_node_mapper.find(subtree_id);
        if (sf_node_it != _subtree_to_sf_node_mapper.end()) {
            return std::optional<NodeId>(sf_node_it->second);
        } else {
            return std::optional<NodeId>{};
        }
    }

private:
    TsNodeToSubtreeMapper const& _ts_node_to_subtree_mapper;
    SubtreeToSfNodeMapper const& _subtree_to_sf_node_mapper;
};
} // namespace sfkit::graph
