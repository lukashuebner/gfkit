#pragma once

#include "tskit/trees.h"
#include <string>

#include "sfkit/graph/primitives.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/samples/primitives.hpp"
#include "sfkit/sequence/Mutation.hpp"
#include "sfkit/sequence/Sequence.hpp"
#include "sfkit/tskit/MutationsView.hpp"

namespace sfkit::tskit {

using ::sfkit::graph::EdgeId;
using ::sfkit::graph::NodeId;
using ::sfkit::graph::TreeId;
using ::sfkit::samples::SampleId;
using ::sfkit::samples::SampleSet;
using ::sfkit::sequence::MutationId;
using ::sfkit::sequence::SiteId;

// TODO Rename to sfkit::tskit::TreeSequence
class TSKitTreeSequence {
public:
    TSKitTreeSequence(std::string const& trees_file);
    TSKitTreeSequence(tsk_treeseq_t const& tree_sequence);
    ~TSKitTreeSequence();

    TSKitTreeSequence(TSKitTreeSequence const& other)            = delete;
    TSKitTreeSequence& operator=(TSKitTreeSequence const& other) = delete;
    TSKitTreeSequence(TSKitTreeSequence&& other) noexcept;

    TSKitTreeSequence& operator=(TSKitTreeSequence&& other);

    [[nodiscard]] NodeId      num_nodes() const;
    [[nodiscard]] EdgeId      num_edges() const;
    [[nodiscard]] SiteId      num_sites() const;
    [[nodiscard]] MutationId  num_mutations() const;
    [[nodiscard]] std::size_t num_populations() const;
    [[nodiscard]] std::size_t num_individuals() const;
    [[nodiscard]] TreeId      num_trees() const;
    [[nodiscard]] SampleId    num_samples() const;
    [[nodiscard]] double      sequence_length() const;

    tsk_treeseq_t&       underlying();
    tsk_treeseq_t const& underlying() const;
    char const*          file_uuid() const;

    [[nodiscard]] bool sample_ids_are_consecutive() const;
    [[nodiscard]] bool is_sample(tsk_id_t u) const;
    [[nodiscard]] bool is_discrete_genome() const;

    [[nodiscard]] double diversity() const;
    [[nodiscard]] double diversity(SampleSet const& samples) const;

    [[nodiscard]] double divergence(SampleSet const& sample_set_1, SampleSet const& sample_set_2) const;

    [[nodiscard]] double f2(SampleSet const& sample_set_1, SampleSet const& sample_set_2) const;
    [[nodiscard]] double
    f3(SampleSet const& sample_set_1, SampleSet const& sample_set_2, SampleSet const& sample_set_3) const;
    [[nodiscard]] double
    f4(SampleSet const& sample_set_1,
       SampleSet const& sample_set_2,
       SampleSet const& sample_set_3,
       SampleSet const& sample_set_4) const;

    [[nodiscard]] double              num_segregating_sites() const;
    [[nodiscard]] std::vector<double> allele_frequency_spectrum() const;

    [[nodiscard]] TskMutationView               mutations() const;
    [[nodiscard]] std::span<tsk_site_t const>   sites() const;
    [[nodiscard]] std::span<double const> const breakpoints() const;

    [[nodiscard]] double position_of(tsk_id_t site_id) const;

private:
    std::string   _trees_file;
    tsk_treeseq_t _tree_sequence;

    bool tskit_noerr(int ret) const;
};

} // namespace sfkit::tskit
