#pragma once

#include <array>
#include <vector>

#include <kassert/kassert.hpp>

#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/samples/sample-set.hpp"
// TODO Unify naming of files
#include "tdt/assertion_levels.hpp"
#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/samples/sample-set.hpp"

template <typename T>
concept NumSamplesBelowAccessorT = requires(T t) {
    { t.num_samples_below(NodeId{}) } -> std::convertible_to<SampleId>;
    { t(NodeId{}) } -> std::convertible_to<SampleId>;
    { t[NodeId{}] } -> std::convertible_to<SampleId>;
    { t.num_nodes_in_dag() } -> std::convertible_to<SampleId>;
    { t.num_samples_in_dag() } -> std::convertible_to<SampleId>;
    { t.num_samples_in_sample_set() } -> std::convertible_to<SampleId>;
};

class NumSamplesBelow {
public:
    NumSamplesBelow(EdgeListGraph const& dag_postorder_edges, SampleSet const& samples)
        : _dag(dag_postorder_edges),
          _num_samples_in_sample_set(samples.popcount()) {
        compute(samples);
    }

    [[nodiscard]] SampleId num_samples_below(NodeId node_id) const {
        KASSERT(node_id < _subtree_sizes.size(), "Subtree ID out of bounds.", tdt::assert::light);
        return _subtree_sizes[node_id];
    }

    [[nodiscard]] SampleId operator()(NodeId node_id) const {
        return this->num_samples_below(node_id);
    }

    [[nodiscard]] SampleId operator[](NodeId node_id) const {
        return this->num_samples_below(node_id);
    }

    [[nodiscard]] SampleId num_nodes_in_dag() const {
        return asserting_cast<SampleId>(_dag.num_nodes());
    }

    [[nodiscard]] SampleId num_samples_in_dag() const {
        return asserting_cast<SampleId>(_dag.num_leaves());
    }

    [[nodiscard]] SampleId num_samples_in_sample_set() const {
        return _num_samples_in_sample_set;
    }

private:
    EdgeListGraph const&  _dag; // As a post-order sorted edge list
    SampleId              _num_samples_in_sample_set;
    std::vector<SampleId> _subtree_sizes;

    // void compute(SampleSet const& samples) {
    //     KASSERT(_dag.check_postorder(), "DAG edges are not post-ordered.", tdt::assert::normal);
    //     KASSERT(
    //         _dag.num_leaves() <= samples.overall_num_samples(),
    //         "Number of leaves in the DAG is greater than he number of overall samples representable in the subtree "
    //         "size object. (NOT the number of samples actually in the SampleSet)",
    //         tdt::assert::light
    //     );

    //     KASSERT(_dag.num_nodes() >= _dag.num_leaves(), "DAG has less nodes than leaves.", tdt::assert::light);
    //     KASSERT(_subtree_sizes.size() == 0ul, "Subtree sizes already computed.", tdt::assert::light);
    //     _subtree_sizes.resize(_dag.num_nodes(), 0);
    //     for (auto&& sample: samples) {
    //         _subtree_sizes[sample] = 1;
    //     }

    //     for (auto& edge: _dag) {
    //         _subtree_sizes[edge.from()] += _subtree_sizes[edge.to()];
    //         KASSERT(
    //             _subtree_sizes[edge.from()] <= samples.popcount(),
    //             "Number of samples below a node exceeds the number of samples in the tree sequence.",
    //             tdt::assert::light
    //         );
    //     }
    // }

    void compute(SampleSet const& samples) {
        KASSERT(_dag.check_postorder(), "DAG edges are not post-ordered.", tdt::assert::normal);
        KASSERT(
            _dag.num_leaves() <= samples.overall_num_samples(),
            "Number of leaves in the DAG is greater than he number of overall samples representable in the subtree "
            "size object. (NOT the number of samples actually in the SampleSet)",
            tdt::assert::light
        );

        KASSERT(_dag.num_nodes() >= _dag.num_leaves(), "DAG has less nodes than leaves.", tdt::assert::light);
        KASSERT(_subtree_sizes.size() == 0ul, "Subtree sizes already computed.", tdt::assert::light);
        _subtree_sizes.resize(_dag.num_nodes(), 0);
        for (auto&& sample: samples) {
            _subtree_sizes[sample] = 1;
        }

        constexpr auto num_edges_to_prefetch = 128;
        auto           prefetch_it           = _dag.begin();
        for (size_t i = 0; i < num_edges_to_prefetch && prefetch_it != _dag.end(); i++) {
            __builtin_prefetch(&_subtree_sizes[prefetch_it->from()], 0, 3);
            __builtin_prefetch(&_subtree_sizes[prefetch_it->to()], 0, 3);
            prefetch_it++;
        }

        auto work_it = _dag.begin();
        while (prefetch_it != _dag.end()) {
            _subtree_sizes[work_it->from()] += _subtree_sizes[work_it->to()];
            KASSERT(
                _subtree_sizes[work_it->from()] <= samples.popcount(),
                "Number of samples below a node exceeds the number of samples in the tree sequence.",
                tdt::assert::light
            );
            work_it++;
            __builtin_prefetch(&_subtree_sizes[prefetch_it->from()], 0, 3);
            __builtin_prefetch(&_subtree_sizes[prefetch_it->to()], 0, 3);
            prefetch_it++;
        }

        while (work_it != _dag.end()) {
            _subtree_sizes[work_it->from()] += _subtree_sizes[work_it->to()];
            KASSERT(
                _subtree_sizes[work_it->from()] <= samples.popcount(),
                "Number of samples below a node exceeds the number of samples in the tree sequence.",
                tdt::assert::light
            );
            work_it++;
        }
    }
};

class NumSamplesBelowTwo {
public:
    NumSamplesBelowTwo(EdgeListGraph const& dag_postorder_edges, SampleSet const& samples_1, SampleSet const& samples_2)
        : _dag(dag_postorder_edges),
          _num_samples_in_sample_set_1(samples_1.popcount()),
          _num_samples_in_sample_set_2(samples_2.popcount()) {
        compute(samples_1, samples_2);
    }

    [[nodiscard]] SampleId num_samples_below(NodeId node_id, uint8_t sample_set) const {
        KASSERT(node_id < _subtree_sizes.size(), "Subtree ID out of bounds.", tdt::assert::light);
        KASSERT(sample_set == 0 || sample_set == 1, "Sample set ID invalid.", tdt::assert::light);

        if (sample_set == 0) {
            return _subtree_sizes[node_id] & _sample_set_1_mask;
        } else {
            return asserting_cast<SampleId>(_subtree_sizes[node_id] >> _sample_set_2_shift);
        }
    }

    [[nodiscard]] SampleId operator()(NodeId node_id, uint8_t sample_set_id) const {
        return this->num_samples_below(node_id, sample_set_id);
    }

    [[nodiscard]] SampleId num_nodes_in_dag() const {
        return asserting_cast<SampleId>(_dag.num_nodes());
    }

    [[nodiscard]] SampleId num_samples_in_dag() const {
        return asserting_cast<SampleId>(_dag.num_leaves());
    }

    [[nodiscard]] SampleId num_samples_in_sample_set(uint8_t sample_set_id) const {
        if (sample_set_id == 0) {
            return _num_samples_in_sample_set_1;
        } else {
            return _num_samples_in_sample_set_2;
        }
    }

private:
    static_assert(
        2 * sizeof(SampleId) == sizeof(size_t), "The number of bits in SampleId is not the half of that of size_t."
    );

    EdgeListGraph const& _dag; // As a post-order sorted edge list
    SampleId             _num_samples_in_sample_set_1;
    SampleId             _num_samples_in_sample_set_2;
    std::vector<size_t>  _subtree_sizes;

    static constexpr size_t _sample_set_1_mask  = 0x00000000FFFFFFFF;
    static constexpr size_t _sample_set_2_mask  = 0xFFFFFFFF00000000;
    static constexpr size_t _sample_set_1_shift = 0;
    static constexpr size_t _sample_set_2_shift = 32;

    void compute(SampleSet const& samples_1, SampleSet const& samples_2) {
        KASSERT(_dag.check_postorder(), "DAG edges are not post-ordered.", tdt::assert::normal);
        KASSERT(
            _dag.num_leaves() <= samples_1.overall_num_samples(),
            "Number of leaves in the DAG is greater than he number of overall samples representable in the subtree "
            "size object. (NOT the number of samples actually in the SampleSet)",
            tdt::assert::light
        );
        KASSERT(
            samples_1.overall_num_samples() == samples_2.overall_num_samples(),
            "Sample sets have different number of overall representable samples. (NOT the number of samples "
            "actually "
            "in the SampleSet)",
            tdt::assert::light
        );

        KASSERT(_dag.num_nodes() >= _dag.num_leaves(), "DAG has less nodes than leaves.", tdt::assert::light);
        KASSERT(_subtree_sizes.size() == 0ul, "Subtree sizes already computed.", tdt::assert::light);
        _subtree_sizes.resize(_dag.num_nodes(), 0);
        for (auto&& sample: samples_1) {
            _subtree_sizes[sample] = 1;
        }
        for (auto&& sample: samples_2) {
            _subtree_sizes[sample] |= 1ul << _sample_set_2_shift;
        }

        constexpr auto num_edges_to_prefetch = 128;
        auto           prefetch_it           = _dag.begin();
        for (size_t i = 0; i < num_edges_to_prefetch && prefetch_it != _dag.end(); i++) {
            __builtin_prefetch(&_subtree_sizes[prefetch_it->from()], 0, 3);
            __builtin_prefetch(&_subtree_sizes[prefetch_it->to()], 0, 3);
            prefetch_it++;
        }

        auto work_it = _dag.begin();
        while (prefetch_it != _dag.end()) {
            // Add both counts bit-parallel
            _subtree_sizes[work_it->from()] += _subtree_sizes[work_it->to()];
            KASSERT(
                (_subtree_sizes[work_it->from()] & _sample_set_1_mask) <= samples_1.popcount(),
                "Number of samples below a node exceeds the number of samples in the tree sequence.",
                tdt::assert::light
            );
            KASSERT(
                (_subtree_sizes[work_it->from()] >> _sample_set_2_shift) <= samples_1.popcount(),
                "Number of samples below a node exceeds the number of samples in the tree sequence.",
                tdt::assert::light
            );
            work_it++;
            __builtin_prefetch(&_subtree_sizes[prefetch_it->from()], 0, 3);
            __builtin_prefetch(&_subtree_sizes[prefetch_it->to()], 0, 3);
            prefetch_it++;
        }

        while (work_it != _dag.end()) {
            _subtree_sizes[work_it->from()] += _subtree_sizes[work_it->to()];
            KASSERT(
                (_subtree_sizes[work_it->from()] & _sample_set_1_mask) <= samples_1.popcount(),
                "Number of samples below a node exceeds the number of samples in the tree sequence.",
                tdt::assert::light
            );
            KASSERT(
                (_subtree_sizes[work_it->from()] >> _sample_set_2_shift) <= samples_1.popcount(),
                "Number of samples below a node exceeds the number of samples in the tree sequence.",
                tdt::assert::light
            );
            work_it++;
        }
    }
};

class NumSamplesBelowTwoAccessor {
public:
    NumSamplesBelowTwoAccessor(std::shared_ptr<NumSamplesBelowTwo> const& num_samples_below, uint8_t sample_set_id)
        : _num_samples_below(num_samples_below),
          _sample_set_id(sample_set_id) {}

    [[nodiscard]] SampleId num_samples_below(NodeId node_id) const {
        return _num_samples_below->operator()(node_id, _sample_set_id);
    }

    [[nodiscard]] SampleId operator()(NodeId node_id) const {
        return this->num_samples_below(node_id);
    }

    [[nodiscard]] SampleId operator[](NodeId node_id) const {
        return this->num_samples_below(node_id);
    }

    [[nodiscard]] SampleId num_nodes_in_dag() const {
        return _num_samples_below->num_nodes_in_dag();
    }

    [[nodiscard]] SampleId num_samples_in_dag() const {
        return _num_samples_below->num_samples_in_dag();
    }

    [[nodiscard]] SampleId num_samples_in_sample_set() const {
        return _num_samples_below->num_samples_in_sample_set(_sample_set_id);
    }

private:
    std::shared_ptr<NumSamplesBelowTwo> _num_samples_below;
    uint8_t                             _sample_set_id;
};

// TODO Use template magic to build on class which can handle all bit-widths and number of sample sets
class NumSamplesBelowFour {
public:
    NumSamplesBelowFour(
        EdgeListGraph const& dag_postorder_edges,
        SampleSet const&     samples_1,
        SampleSet const&     samples_2,
        SampleSet const&     samples_3,
        SampleSet const&     samples_4
    )
        : _dag(dag_postorder_edges),
          _num_samples_in_sample_set(
              {samples_1.popcount(), samples_2.popcount(), samples_3.popcount(), samples_4.popcount()}
          ) {
        for (auto&& num_samples: _num_samples_in_sample_set) {
            KASSERT(
                num_samples <= samples_1.overall_num_samples(),
                "Number of samples in sample set is greater than he number of overall samples representable in the "
                "subtree size object. (NOT the number of samples actually in the SampleSet)",
                tdt::assert::light
            );
        }
        compute(samples_1, samples_2, samples_3, samples_4);
    }

    [[nodiscard]] SampleId num_samples_below(NodeId node_id, uint8_t sample_set) const {
        KASSERT(node_id < _subtree_sizes.size(), "Subtree ID out of bounds.", tdt::assert::light);
        KASSERT(sample_set >= 0 && sample_set <= 3, "Sample set ID invalid.", tdt::assert::light);

        return asserting_cast<SampleId>(
            (_subtree_sizes[node_id] & _sample_set_masks[sample_set]) >> _sample_set_shifts[sample_set]
        );
    }

    [[nodiscard]] SampleId operator()(NodeId node_id, uint8_t sample_set_id) const {
        return this->num_samples_below(node_id, sample_set_id);
    }

    [[nodiscard]] SampleId num_nodes_in_dag() const {
        return asserting_cast<SampleId>(_dag.num_nodes());
    }

    [[nodiscard]] SampleId num_samples_in_dag() const {
        return asserting_cast<SampleId>(_dag.num_leaves());
    }

    [[nodiscard]] SampleId num_samples_in_sample_set(uint8_t sample_set_id) const {
        KASSERT(sample_set_id >= 0 && sample_set_id <= 3, "Sample set ID invalid.", tdt::assert::light);
        return _num_samples_in_sample_set[sample_set_id];
    }

private:
    // TODO Try using uint128_t
    static_assert(
        2 * sizeof(SampleId) == sizeof(size_t), "The number of bits in SampleId is not the half of that of size_t."
    );

    EdgeListGraph const&    _dag; // As a post-order sorted edge list
    std::array<SampleId, 4> _num_samples_in_sample_set;
    std::vector<size_t>     _subtree_sizes;

    static constexpr std::array<size_t, 4> _sample_set_shifts = {0, 16, 32, 48};
    static constexpr std::array<size_t, 4> _sample_set_masks  = {
         0x000000000000FFFF, 0x0000000FFFF0000, 0x0000FFFF00000000, 0xFFFF000000000000};

    void compute(
        SampleSet const& samples_1, SampleSet const& samples_2, SampleSet const& samples_3, SampleSet const& samples_4
    ) {
        KASSERT(_dag.check_postorder(), "DAG edges are not post-ordered.", tdt::assert::normal);
        KASSERT(
            _dag.num_leaves() <= samples_1.overall_num_samples(),
            "Number of leaves in the DAG is greater than he number of overall samples representable in the subtree "
            "size object. (NOT the number of samples actually in the SampleSet)",
            tdt::assert::light
        );
        KASSERT(
            samples_1.overall_num_samples() == samples_2.overall_num_samples(),
            "Sample sets have different number of overall representable samples. (NOT the number of samples "
            "actually "
            "in the SampleSet)",
            tdt::assert::light
        );

        KASSERT(_dag.num_nodes() >= _dag.num_leaves(), "DAG has less nodes than leaves.", tdt::assert::light);
        KASSERT(_subtree_sizes.size() == 0ul, "Subtree sizes already computed.", tdt::assert::light);
        _subtree_sizes.resize(_dag.num_nodes(), 0);

        for (auto&& sample: samples_1) {
            _subtree_sizes[sample] = 1ul;
        }
        for (auto&& sample: samples_2) {
            _subtree_sizes[sample] |= 1ul << _sample_set_shifts[1];
        }
        for (auto&& sample: samples_3) {
            _subtree_sizes[sample] |= 1ul << _sample_set_shifts[2];
        }
        // TODO Use 1-based indexing for sample sets
        for (auto&& sample: samples_4) {
            _subtree_sizes[sample] |= 1ul << _sample_set_shifts[3];
        }

        constexpr auto num_edges_to_prefetch = 128;
        auto           prefetch_it           = _dag.begin();
        for (size_t i = 0; i < num_edges_to_prefetch && prefetch_it != _dag.end(); i++) {
            __builtin_prefetch(&_subtree_sizes[prefetch_it->from()], 0, 3);
            __builtin_prefetch(&_subtree_sizes[prefetch_it->to()], 0, 3);
            prefetch_it++;
        }

        auto work_it = _dag.begin();
        while (prefetch_it != _dag.end()) {
            // Add all four bit counts bit-parallel
            _subtree_sizes[work_it->from()] += _subtree_sizes[work_it->to()];
            for (auto mask: _sample_set_masks) {
                KASSERT(
                    (_subtree_sizes[work_it->from()] & mask) <= samples_1.popcount(),
                    "Number of samples below a node exceeds the number of samples in the tree sequence.",
                    tdt::assert::light
                );
            }
            work_it++;
            __builtin_prefetch(&_subtree_sizes[prefetch_it->from()], 0, 3);
            __builtin_prefetch(&_subtree_sizes[prefetch_it->to()], 0, 3);
            prefetch_it++;
        }

        while (work_it != _dag.end()) {
            _subtree_sizes[work_it->from()] += _subtree_sizes[work_it->to()];
            for (auto mask: _sample_set_masks) {
                KASSERT(
                    (_subtree_sizes[work_it->from()] & mask) <= samples_1.popcount(),
                    "Number of samples below a node exceeds the number of samples in the tree sequence.",
                    tdt::assert::light
                );
            }
            work_it++;
        }
    }
};

// TODO Use a single Accessor for One, Two, and Four
class NumSamplesBelowFourAccessor {
public:
    NumSamplesBelowFourAccessor(std::shared_ptr<NumSamplesBelowFour> const& num_samples_below, uint8_t sample_set_id)
        : _num_samples_below(num_samples_below),
          _sample_set_id(sample_set_id) {}

    [[nodiscard]] SampleId num_samples_below(NodeId node_id) const {
        return _num_samples_below->operator()(node_id, _sample_set_id);
    }

    [[nodiscard]] SampleId operator()(NodeId node_id) const {
        return this->num_samples_below(node_id);
    }

    [[nodiscard]] SampleId operator[](NodeId node_id) const {
        return this->num_samples_below(node_id);
    }

    [[nodiscard]] SampleId num_nodes_in_dag() const {
        return _num_samples_below->num_nodes_in_dag();
    }

    [[nodiscard]] SampleId num_samples_in_dag() const {
        return _num_samples_below->num_samples_in_dag();
    }

    [[nodiscard]] SampleId num_samples_in_sample_set() const {
        return _num_samples_below->num_samples_in_sample_set(_sample_set_id);
    }

private:
    std::shared_ptr<NumSamplesBelowFour> _num_samples_below;
    uint8_t                              _sample_set_id;
};
