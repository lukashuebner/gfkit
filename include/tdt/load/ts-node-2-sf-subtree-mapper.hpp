#pragma once

#include <algorithm>
#include <vector>

#include <kassert/kassert.hpp>
#include <tskit/core.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>
#pragma GCC diagnostic pop

#include "tdt/graph/common.hpp"
#include "tdt/utils/concepts.hpp"

// class TsNode2SfSubtreeMapper {
// public:
//     void reserve(size_t size) {
//         _map.reserve(size);
//     }
//
//     void reset() {
//         _map.clear();
//     }
//
//     NodeId& operator[](tsk_id_t const ts_node_id) {
//         return map(ts_node_id);
//     }
//
//     NodeId& map(tsk_id_t const ts_node_id) {
//         return _map[ts_node_id];
//     }
//
//     NodeId& at(tsk_id_t const ts_node_id) {
//         return _map.at(ts_node_id);
//     }
//
//     NodeId const& at(tsk_id_t const ts_node_id) const {
//         return _map.at(ts_node_id);
//     }
//
//     template <IterableInput IterableInput>
//     std::vector<NodeId> map(IterableInput input) const {
//         std::vector<NodeId> result;
//         result.reserve(input.size());
//         std::transform(input.begin(), input.end(), std::back_inserter(result), [this](auto&& tsk_node) {
//             return map(tsk_node);
//         });
//         return result;
//     }
//
// private:
//     mutable tsl::hopscotch_map<tsk_id_t, NodeId> _map;
// };

// class TsToSfNodeMapper {
// public:
//     void reserve(size_t size) {
//         _map.reserve(size);
//         _requested.reserve(size);
//     }

//     void reset() {
//         _map.clear();
//         _requested.clear();
//     }

//     void insert_if_requested(tsk_id_t ts_node_id, NodeId dag_node_id) {
//         if (is_requested(ts_node_id)) [[unlikely]] {
//             insert(ts_node_id, dag_node_id);
//         }
//     }

//     void insert(tsk_id_t ts_node_id, NodeId node_id) {
//         KASSERT(is_requested(ts_node_id), "The vertex " << ts_node_id << "was not requested", tdt::assert::normal);
//         _map.insert({ts_node_id, node_id});
//     }

//     NodeId const& at(tsk_id_t const ts_node_id) const {
//         KASSERT(is_requested(ts_node_id), "The vertex " << ts_node_id << "was not requested", tdt::assert::normal);
//         return _map.at(ts_node_id);
//     }

//     void request(tsk_id_t const ts_node_id) {
//         _requested.insert(ts_node_id);
//     }

//     bool is_requested(tsk_id_t const ts_node_id) const {
//         return _requested.contains(ts_node_id);
//     }

//     std::size_t size() const {
//         return _map.size();
//     }

//     std::size_t num_requested() const {
//         return _requested.size();
//     }

// private:
//     mutable tsl::hopscotch_map<tsk_id_t, NodeId> _map;
//     tsl::hopscotch_set<tsk_id_t>                 _requested;
// };

// TODO Write a general mapper which also supports chaining
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

private:
    TsNodeToSubtreeMapper const& _ts_node_to_subtree_mapper;
    SubtreeToSfNodeMapper const& _subtree_to_sf_node_mapper;
};
