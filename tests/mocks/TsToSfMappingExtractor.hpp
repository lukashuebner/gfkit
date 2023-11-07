#pragma once

#include "tskit/core.h"
#include <unordered_map>
#include <vector>

#include "sfkit/checking_casts.hpp"
#include "sfkit/graph/common.hpp"
#include "sfkit/load/SubtreeHashToNodeMapper.hpp"
#include "sfkit/load/SubtreeHasher.hpp"
#include "sfkit/load/TsToSfNodeMapper.hpp"

class Ts2SfMappingExtractor {
public:
    using SingleTreeMapper = TsToSfNodeMapper<std::vector<SubtreeHash>, SubtreeHashToNodeMapper>;

    Ts2SfMappingExtractor(TreeId num_trees, NodeId num_nodes) : _num_nodes(num_nodes) {
        _mappings.resize(num_trees);
    }

    bool process_mutations(TreeId tree_id, SingleTreeMapper const& mapper) {
        _process_mutations_callcnt++;

        for (tsk_id_t ts_node = 0; ts_node < asserting_cast<tsk_id_t>(_num_nodes); ts_node++) {
            auto const dag_node = mapper.try_map(asserting_cast<size_t>(ts_node));
            if (dag_node) {
                _mappings[tree_id][ts_node] = *dag_node;
            }
        }

        return true;
    }

    void finalize() {
        _finalize_called = true;
    }

    NodeId operator()(TreeId const tree, tsk_id_t const ts_id) {
        KASSERT(tree < _mappings.size(), "Tree ID out of bounds.", sfkit::assert::light);
        KASSERT(_mappings[tree].find(ts_id) != _mappings[tree].end(), "TS node ID not found.", sfkit::assert::normal);
        return _mappings[tree][ts_id];
    }

    template <IterableInput IterableInput>
    std::vector<NodeId> operator()(TreeId const tree, IterableInput input) {
        std::vector<NodeId> result;
        result.reserve(input.size());
        std::transform(input.begin(), input.end(), std::back_inserter(result), [this, tree](tsk_id_t tsk_node) {
            return operator()(tree, tsk_node);
        });
        return result;
    }

    bool finalize_called() {
        return _finalize_called;
    }

    size_t process_mutations_callcnt() {
        return _process_mutations_callcnt;
    }

private:
    using AllTreesMap = std::vector<std::unordered_map<tsk_id_t, NodeId>>;

    bool        _finalize_called           = false;
    size_t      _process_mutations_callcnt = 0;
    NodeId      _num_nodes;
    AllTreesMap _mappings;
};
