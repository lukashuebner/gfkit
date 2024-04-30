#pragma once

#include <string>
#include <vector>

#include <kassert/kassert.hpp>
#include <tskit.h>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/bp/BPCompressedForest.hpp"
#include "sfkit/bp/BPForestCompressor.hpp"
#include "sfkit/dag/DAGCompressedForest.hpp"
#include "sfkit/dag/DAGForestCompressor.hpp"
#include "sfkit/graph/ForestCompressor.hpp"
#include "sfkit/samples/NumSamplesBelowFactory.hpp"
#include "sfkit/sequence/AlleleFrequencies.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"
#include "sfkit/sequence/GenomicSequenceFactory.hpp"
#include "sfkit/stats/AlleleFrequencySpectrum.hpp"
#include "sfkit/stats/Divergence.hpp"
#include "sfkit/stats/Diversity.hpp"
#include "sfkit/stats/Fst.hpp"
#include "sfkit/stats/LCA.hpp"
#include "sfkit/stats/NumSegregatingSites.hpp"
#include "sfkit/stats/PattersonsF.hpp"
#include "sfkit/stats/TajimasD.hpp"
#include "sfkit/tskit/tskit.hpp"
#include "sfkit/utils/always_false_v.hpp"
#include "sfkit/utils/tuple_transform.hpp"

namespace sfkit {

using sfkit::bp::BPCompressedForest;
using sfkit::dag::DAGCompressedForest;
using sfkit::graph::ForestCompressor;
using sfkit::sequence::AlleleFrequencies;
using sfkit::tskit::TSKitTreeSequence;
using namespace sfkit::samples;
using namespace sfkit::sequence;

// TODO Re-think the interface for this (esp. construction)
template <typename CompressedForest, typename PerfectAllelicStateHasher>
class SuccinctForest {
public:
    SuccinctForest(CompressedForest&& forest, GenomicSequence&& genomic_sequence)
        : _forest(forest),
          _sequence(genomic_sequence) {}

    explicit SuccinctForest(tsk_treeseq_t& ts_tree_sequence) {
        TSKitTreeSequence tree_sequence(ts_tree_sequence);
        _init(tree_sequence);
    }

    explicit SuccinctForest(TSKitTreeSequence& tree_sequence) {
        _init(tree_sequence);
    }

    explicit SuccinctForest(std::string const& filename) {
        TSKitTreeSequence tree_sequence(filename);
        _init(tree_sequence);
    }

    // TODO Pass by reference
    auto allele_frequencies(SampleSet const samples) {
        using NumSamplesBelowT = NumSamplesBelow<CompressedForest, 1, SampleId>;
        auto num_samples_below = NumSamplesBelowAccessor<NumSamplesBelowT>(
            std::make_shared<NumSamplesBelowT>(
                _forest,
                std::array<std::reference_wrapper<SampleSet const>, 1>{std::cref(samples)}
            ),
            0
        );
        return allele_frequencies(num_samples_below);
    }

    template <NumSamplesBelowAccessorC NumSamplesBelowAccessorT>
    auto allele_frequencies(NumSamplesBelowAccessorT const& num_samples_below) {
        return AlleleFrequencies<CompressedForest, PerfectAllelicStateHasher, NumSamplesBelowAccessorT>(
            _forest,
            _sequence,
            num_samples_below
        );
    }

    // TODO Pass by reference
    auto allele_frequencies(SampleSet const samples_0, SampleSet const samples_1) {
        return sfkit::utils::tuple_transform(
            [this](auto&& e) { return allele_frequencies(e); },
            NumSamplesBelowFactory::build(_forest, samples_0, samples_1)
        );
    }

    // TODO Pass by reference
    template <typename NumSamplesBelowBaseType = SampleId>
    auto allele_frequencies(SampleSet const samples_0, SampleSet const samples_1, SampleSet const samples_2) {
        SampleSet empty_sample_set = SampleSet(samples_0.overall_num_samples());
        auto [num_samples_below_0, num_samples_below_1, num_samples_below_2, dummy] =
            NumSamplesBelowFactory::build<CompressedForest, NumSamplesBelowBaseType>(
                _forest,
                samples_0,
                samples_1,
                samples_2,
                empty_sample_set
            );
        return std::tuple(
            allele_frequencies(num_samples_below_0),
            allele_frequencies(num_samples_below_1),
            allele_frequencies(num_samples_below_2)
        );
    }

    // TODO Pass by reference
    template <typename NumSamplesBelowBaseType = SampleId>
    auto allele_frequencies(
        SampleSet const samples_0, SampleSet const samples_1, SampleSet const samples_2, SampleSet const samples_3
    ) {
        return sfkit::utils::tuple_transform(
            [this](auto&& e) { return allele_frequencies(e); },
            NumSamplesBelowFactory::build<CompressedForest, NumSamplesBelowBaseType>(
                _forest,
                samples_0,
                samples_1,
                samples_2,
                samples_3
            )
        );
    }

    // TODO Make this const
    [[nodiscard]] double diversity() {
        return diversity(_forest.all_samples());
    }

    template <typename AlleleFrequencies>
    [[nodiscard]] double diversity(SampleId num_samples, AlleleFrequencies allele_freqs) {
        return stats::Diversity::diversity(num_samples, allele_freqs);
    }

    // TODO Pass by reference
    [[nodiscard]] double diversity(SampleSet const sample_set) {
        SampleId const num_samples = sample_set.popcount();
        auto const     freqs       = allele_frequencies(sample_set);
        return diversity(num_samples, freqs);
    }

    [[nodiscard]] auto allele_frequency_spectrum() {
        return allele_frequency_spectrum(_forest.all_samples());
    }

    // TODO Pass by reference
    [[nodiscard]] auto allele_frequency_spectrum(SampleSet sample_set) {
        return sfkit::stats::AlleleFrequencySpectrum(allele_frequencies(sample_set));
    }

    template <typename AlleleFrequenciesT>
    [[nodiscard]] double divergence(
        SampleId           num_samples_0,
        AlleleFrequenciesT allele_frequencies_0,
        SampleId           num_samples_1,
        AlleleFrequenciesT allele_frequencies_1
    ) {
        return stats::Divergence::divergence(num_samples_0, allele_frequencies_0, num_samples_1, allele_frequencies_1);
    }

    // TODO Pass by reference?
    [[nodiscard]] double divergence(SampleSet const sample_set_0, SampleSet const sample_set_1) {
        auto [allele_freqs_0, allele_freqs_1] = allele_frequencies(sample_set_0, sample_set_1);
        auto const num_samples_0              = sample_set_0.popcount();
        auto const num_samples_1              = sample_set_1.popcount();
        return divergence(num_samples_0, allele_freqs_0, num_samples_1, allele_freqs_1);
    }

    // TODO Pass by reference?
    [[nodiscard]] double f2(SampleSet const sample_set_0, SampleSet const sample_set_1) {
        auto const [allele_freqs_0, allele_freqs_1] = allele_frequencies(sample_set_0, sample_set_1);
        return stats::PattersonsF::f2(allele_freqs_0, allele_freqs_1);
    }

    // TODO Pass by reference?
    [[nodiscard]] double f3(SampleSet const samples_0, SampleSet const samples_1, SampleSet const samples_2) {
        if (samples_0.popcount() <= UINT16_MAX && samples_1.popcount() <= UINT16_MAX
            && samples_2.popcount() <= UINT16_MAX) [[likely]] {
            auto const [allele_freqs_0, allele_freqs_1, allele_freqs_2] =
                allele_frequencies<uint16_t>(samples_0, samples_1, samples_2);
            return stats::PattersonsF::f3(allele_freqs_0, allele_freqs_1, allele_freqs_2);
        } else {
            auto const [allele_freqs_0, allele_freqs_1, allele_freqs_2] =
                allele_frequencies<SampleId>(samples_0, samples_1, samples_2);
            return stats::PattersonsF::f3(allele_freqs_0, allele_freqs_1, allele_freqs_2);
        }
    }

    // TODO Pass by reference?
    [[nodiscard]] double
    f4(SampleSet const samples_0, SampleSet const samples_1, SampleSet const samples_2, SampleSet const samples_3) {
        if (samples_0.popcount() <= UINT16_MAX && samples_1.popcount() <= UINT16_MAX
            && samples_2.popcount() <= UINT16_MAX && samples_3.popcount() <= UINT16_MAX) [[likely]] {
            auto const [allele_freqs_0, allele_freqs_1, allele_freqs_2, allele_freqs_3] =
                allele_frequencies<uint16_t>(samples_0, samples_1, samples_2, samples_3);
            return stats::PattersonsF::f4(allele_freqs_0, allele_freqs_1, allele_freqs_2, allele_freqs_3);
        } else {
            auto const [allele_freqs_0, allele_freqs_1, allele_freqs_2, allele_freqs_3] =
                allele_frequencies<SampleId>(samples_0, samples_1, samples_2, samples_3);
            return stats::PattersonsF::f4(allele_freqs_0, allele_freqs_1, allele_freqs_2, allele_freqs_3);
        }
    }

    [[nodiscard]] std::vector<NodeId> lca(SampleId const u, SampleId const v) {
        SampleSet samples(_forest.num_samples());
        samples.add(u);
        samples.add(v);

        return lca(samples);
    }

    [[nodiscard]] std::vector<NodeId> lca(SampleSet const& samples) {
        if constexpr (std::is_same_v<CompressedForest, DAGCompressedForest>) {
            stats::DAGLowestCommonAncestor lca(_forest.postorder_edges());
            return lca.lca(samples);
        } else {
            KASSERT(false, "LCA not implemented for BP variant of SuccinctForest", sfkit::assert::light);

            return std::vector<NodeId>(num_trees(), graph::INVALID_NODE_ID);
        }
    }

    // TODO Make this const
    template <typename AlleleFrequencies>
    [[nodiscard]] SiteId num_segregating_sites(SampleId num_samples, AlleleFrequencies allele_frequencies) {
        return sfkit::stats::NumSegregatingSites::num_segregating_sites(num_samples, allele_frequencies);
    }

    // TODO Pass by reference?
    [[nodiscard]] SiteId num_segregating_sites(SampleSet const sample_set) {
        auto const num_samples = sample_set.popcount();
        auto const freqs       = allele_frequencies(sample_set);
        return num_segregating_sites(num_samples, freqs);
    }

    [[nodiscard]] SiteId num_segregating_sites() {
        return num_segregating_sites(_forest.all_samples());
    }

    [[nodiscard]] double tajimas_d() {
        auto const allele_freqs = allele_frequencies(_forest.all_samples());
        return stats::TajimasD::tajimas_d(num_samples(), allele_freqs);
    }

    // This is per sequence length, the other statistics are not
    // TODO Pass by reference?
    [[nodiscard]] double fst(SampleSet const sample_set_0, SampleSet const sample_set_1) {
        auto [allele_freqs_0, allele_freqs_1] = allele_frequencies(sample_set_0, sample_set_1);
        return stats::Fst::fst(_sequence.num_sites(), allele_freqs_0, allele_freqs_1);
    }

    [[nodiscard]] SiteId num_sites() const {
        return _sequence.num_sites();
    }

    [[nodiscard]] SampleId num_samples() const {
        return _forest.num_samples();
    }

    [[nodiscard]] SampleSet all_samples() const {
        return _forest.all_samples();
    }

    [[nodiscard]] GenomicSequence const& sequence() const {
        return _sequence;
    }

    [[nodiscard]] CompressedForest const& forest() const {
        return _forest;
    }

    [[nodiscard]] TreeId num_trees() const {
        return _forest.num_trees();
    }

    [[nodiscard]] MutationId num_mutations() const {
        return _sequence.num_mutations();
    }

    [[nodiscard]] auto num_unique_subtrees() const {
        return _forest.num_unique_subtrees();
    }

    [[nodiscard]] auto num_subtrees_with_mutations() const {
        return _sequence.subtrees_with_mutations().size();
    }

private:
    CompressedForest _forest;
    GenomicSequence  _sequence;

    void _init(TSKitTreeSequence& tree_sequence) {
        ForestCompressor<CompressedForest> forest_compressor(tree_sequence);
        GenomicSequenceFactory             sequence_factory(tree_sequence);
        _forest   = forest_compressor.compress(sequence_factory);
        _sequence = sequence_factory.move_storage();
    }
};

using DAGSuccinctForest        = SuccinctForest<DAGCompressedForest, PerfectDNAHasher>;
using DAGSuccinctForestNumeric = SuccinctForest<DAGCompressedForest, PerfectNumericHasher>;
using BPSuccinctForest         = SuccinctForest<BPCompressedForest, PerfectDNAHasher>;
using BPSuccinctForestNumeric  = SuccinctForest<BPCompressedForest, PerfectNumericHasher>;

} // namespace sfkit
