#pragma once

#include <algorithm>
#include <vector>

#include <tskit/core.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include <tsl/hopscotch_map.h>
#pragma GCC diagnostic pop

#include "tdt/graph/common.hpp"
#include "tdt/utils/concepts.hpp"

class TsNode2SfSubtreeMapper {
public:
    void reserve(size_t size) {
        _map.reserve(size);
    }

    void reset() {
        _map.clear();
    }

    NodeId& operator[](tsk_id_t const ts_node_id) {
        return map(ts_node_id);
    }

    NodeId& map(tsk_id_t const ts_node_id) {
        return _map[ts_node_id];
    }

    NodeId& at(tsk_id_t const ts_node_id) {
        return _map.at(ts_node_id);
    }

    NodeId const& at(tsk_id_t const ts_node_id) const {
        return _map.at(ts_node_id);
    }

    template <IterableInput IterableInput>
    std::vector<NodeId> map(IterableInput input) const {
        std::vector<NodeId> result;
        result.reserve(input.size());
        std::transform(input.begin(), input.end(), std::back_inserter(result), [this](auto&& tsk_node) {
            return map(tsk_node);
        });
        return result;
    }

    std::size_t size() const {
        return _map.size();
    }
private:
    mutable tsl::hopscotch_map<tsk_id_t, NodeId> _map;
};
