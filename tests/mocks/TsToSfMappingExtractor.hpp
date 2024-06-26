#pragma once

#include "tskit/core.h"
#include <unordered_map>
#include <vector>

#include "sfkit/graph/SubtreeHashToNodeMapper.hpp"
#include "sfkit/graph/SubtreeHasher.hpp"
#include "sfkit/graph/TsToSfNodeMapper.hpp"
#include "sfkit/graph/primitives.hpp"
#include "sfkit/utils/checking_casts.hpp"

using namespace sfkit::graph;
using sfkit::utils::asserting_cast;
using sfkit::utils::IterableInput;

class Sf2TsMapper {
public:
    Sf2TsMapper(std::vector<std::unordered_map<tsk_id_t, NodeId>> ts2sf) {
        _sf2ts.resize(ts2sf.size());
        for (TreeId tree_id = 0; tree_id < ts2sf.size(); ++tree_id) {
            for (auto const [key, value]: ts2sf[tree_id]) {
                _sf2ts[tree_id][value] = key;
            }
        }
    }

    tsk_id_t operator()(TreeId const tree, NodeId const sf_id) {
        KASSERT(tree < _sf2ts.size(), "Tree ID out of bounds.", sfkit::assert::light);
        KASSERT(_sf2ts[tree].find(sf_id) != _sf2ts[tree].end(), "TS node ID not found.", sfkit::assert::normal);
        return _sf2ts[tree][sf_id];
    }

private:
    std::vector<std::unordered_map<NodeId, tsk_id_t>> _sf2ts;
};

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

    Sf2TsMapper inverse() {
        return Sf2TsMapper{_mappings};
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
