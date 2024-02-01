#include "sfkit/tskit/TSKitTreeSequence.hpp"

#include <cstddef>
#include <iterator>
#include <numeric>
#include <span>

#include <kassert/kassert.hpp>
#include <tskit.h>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/graph/primitives.hpp"
#include "sfkit/samples/SampleSet.hpp"
#include "sfkit/sequence/Mutation.hpp"
#include "sfkit/sequence/Sequence.hpp"
#include "sfkit/utils/checking_casts.hpp"

namespace sfkit::tskit {

using sfkit::graph::EdgeId;
using sfkit::graph::NodeId;
using sfkit::graph::TreeId;
using sfkit::samples::SampleId;
using sfkit::samples::SampleSet;
using sfkit::sequence::MutationId;
using sfkit::sequence::SiteId;
using sfkit::utils::asserting_cast;

TSKitTreeSequence::TSKitTreeSequence(std::string const& trees_file) : _trees_file(trees_file), owning(true) {
    int ret;

    // Load the tree sequence from the .trees file
    ret = tsk_treeseq_load(&_tree_sequence, _trees_file.c_str(), 0);
    KASSERT(tskit_noerr(ret), "Failed to load tree sequence from the .trees file", sfkit::assert::light);
}

TSKitTreeSequence::TSKitTreeSequence(tsk_treeseq_t const& tree_sequence)
    : _tree_sequence(tree_sequence),
      owning(false) {
}

TSKitTreeSequence::TSKitTreeSequence(tsk_treeseq_t const&& tree_sequence)
    : _tree_sequence(tree_sequence),
      owning(true) {
}

TSKitTreeSequence::~TSKitTreeSequence() {
    if (owning) {
        tsk_treeseq_free(&_tree_sequence);
    }
}

TSKitTreeSequence::TSKitTreeSequence(TSKitTreeSequence&& other) noexcept
    : _tree_sequence(other._tree_sequence),
      owning(other.owning) {
    other._tree_sequence = tsk_treeseq_t();
    other.owning         = false;
}

TSKitTreeSequence& TSKitTreeSequence::operator=(TSKitTreeSequence&& other) {
    if (this != &other) {
        tsk_treeseq_free(&_tree_sequence);
        _tree_sequence       = other._tree_sequence;
        other._tree_sequence = tsk_treeseq_t();
    }
    return *this;
}

bool TSKitTreeSequence::is_owning() const {
    return owning;
}

NodeId TSKitTreeSequence::num_nodes() const {
    return asserting_cast<NodeId>(tsk_treeseq_get_num_nodes(&_tree_sequence));
}

EdgeId TSKitTreeSequence::num_edges() const {
    return asserting_cast<EdgeId>(tsk_treeseq_get_num_edges(&_tree_sequence));
}

SiteId TSKitTreeSequence::num_sites() const {
    return asserting_cast<SiteId>(tsk_treeseq_get_num_sites(&_tree_sequence));
}

MutationId TSKitTreeSequence::num_mutations() const {
    return asserting_cast<MutationId>(tsk_treeseq_get_num_mutations(&_tree_sequence));
}

std::size_t TSKitTreeSequence::num_populations() const {
    return tsk_treeseq_get_num_populations(&_tree_sequence);
}

std::size_t TSKitTreeSequence::num_individuals() const {
    return tsk_treeseq_get_num_individuals(&_tree_sequence);
}

TreeId TSKitTreeSequence::num_trees() const {
    return asserting_cast<TreeId>(tsk_treeseq_get_num_trees(&_tree_sequence));
}

SampleId TSKitTreeSequence::num_samples() const {
    return asserting_cast<SampleId>(tsk_treeseq_get_num_samples(&_tree_sequence));
}

char const* TSKitTreeSequence::file_uuid() const {
    return tsk_treeseq_get_file_uuid(&_tree_sequence);
}

double TSKitTreeSequence::sequence_length() const {
    return tsk_treeseq_get_sequence_length(&_tree_sequence);
}

tsk_treeseq_t& TSKitTreeSequence::underlying() {
    return _tree_sequence;
}

tsk_treeseq_t const& TSKitTreeSequence::underlying() const {
    return _tree_sequence;
}

[[nodiscard]] bool TSKitTreeSequence::sample_ids_are_consecutive() const {
    // According to the tskit documentation: The array is owned by the tree sequence and should not be modified
    // or free’d.
    tsk_id_t const*  samples     = tsk_treeseq_get_samples(&underlying());
    tsk_size_t const num_samples = tsk_treeseq_get_num_samples(&underlying());
    for (tsk_id_t id = 0; id < asserting_cast<tsk_id_t>(num_samples); ++id) {
        if (id != *(samples + id)) {
            return false;
        }
    }
    return true;
}

bool TSKitTreeSequence::is_sample(tsk_id_t u) const {
    return tsk_treeseq_is_sample(&_tree_sequence, u);
}

bool TSKitTreeSequence::is_discrete_genome() const {
    return tsk_treeseq_get_discrete_genome(&_tree_sequence);
}

double TSKitTreeSequence::diversity() const {
    double                pi;
    std::vector<tsk_id_t> samples;
    samples.resize(num_samples());
    std::iota(samples.begin(), samples.end(), 0);
    tsk_size_t sample_set_sizes[] = {num_samples()};

    auto ret = tsk_treeseq_diversity(&_tree_sequence, 1, sample_set_sizes, samples.data(), 0, NULL, TSK_STAT_SITE, &pi);
    KASSERT(ret == 0, "Failed to compute the diversity.", sfkit::assert::light);

    return pi;
}

double TSKitTreeSequence::diversity(SampleSet const& samples) const {
    double pi;

    auto const tsk_samples        = samples.to_tsk_samples();
    tsk_size_t sample_set_sizes[] = {tsk_samples.size()};

    auto ret =
        tsk_treeseq_diversity(&_tree_sequence, 1, sample_set_sizes, tsk_samples.data(), 0, NULL, TSK_STAT_SITE, &pi);
    KASSERT(ret == 0, "Failed to compute the diversity.", sfkit::assert::light);

    return pi;
}

double TSKitTreeSequence::divergence(SampleSet const& sample_set_1, SampleSet const& sample_set_2) const {
    double                divergence;
    constexpr int         num_windows     = 0;
    constexpr int         num_sample_sets = 2;
    std::vector<tsk_id_t> samples_sets;

    for (auto const& sample: sample_set_1) {
        samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
    }
    for (auto const& sample: sample_set_2) {
        samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
    }
    tsk_size_t sample_set_sizes[] = {sample_set_1.popcount(), sample_set_2.popcount()};

    tsk_id_t      set_indexes[]    = {0, 1};
    constexpr int num_index_tuples = 1;

    auto ret = tsk_treeseq_divergence(
        &_tree_sequence,
        num_sample_sets,
        sample_set_sizes,
        samples_sets.data(),
        num_index_tuples,
        set_indexes,
        num_windows,
        NULL,
        TSK_STAT_SITE,
        &divergence
    );
    KASSERT(ret == 0, "Failed to compute the divergence.", sfkit::assert::light);

    return divergence;
}

double TSKitTreeSequence::f4(
    SampleSet const& sample_set_1,
    SampleSet const& sample_set_2,
    SampleSet const& sample_set_3,
    SampleSet const& sample_set_4
) const {
    double                f4;
    constexpr int         num_windows     = 0;
    constexpr int         num_sample_sets = 4;
    std::vector<tsk_id_t> samples_sets;

    for (auto const& sample: sample_set_1) {
        samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
    }
    for (auto const& sample: sample_set_2) {
        samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
    }
    for (auto const& sample: sample_set_3) {
        samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
    }
    for (auto const& sample: sample_set_4) {
        samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
    }
    tsk_size_t sample_set_sizes[] =
        {sample_set_1.popcount(), sample_set_2.popcount(), sample_set_3.popcount(), sample_set_4.popcount()};

    tsk_id_t      set_indexes[]    = {0, 1, 2, 3};
    constexpr int num_index_tuples = 1;

    auto ret = tsk_treeseq_f4(
        &_tree_sequence,
        num_sample_sets,
        sample_set_sizes,
        samples_sets.data(),
        num_index_tuples,
        set_indexes,
        num_windows,
        NULL,
        TSK_STAT_SITE,
        &f4
    );
    KASSERT(ret == 0, "Failed to compute the f4.", sfkit::assert::light);

    return f4;
}

double TSKitTreeSequence::f3(
    SampleSet const& sample_set_1, SampleSet const& sample_set_2, SampleSet const& sample_set_3
) const {
    double                f3;
    constexpr int         num_windows     = 0;
    constexpr int         num_sample_sets = 3;
    std::vector<tsk_id_t> samples_sets;

    for (auto const& sample: sample_set_1) {
        samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
    }
    for (auto const& sample: sample_set_2) {
        samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
    }
    for (auto const& sample: sample_set_3) {
        samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
    }
    tsk_size_t sample_set_sizes[] = {sample_set_1.popcount(), sample_set_2.popcount(), sample_set_3.popcount()};

    tsk_id_t      set_indexes[]    = {0, 1, 2, 3};
    constexpr int num_index_tuples = 1;

    auto ret = tsk_treeseq_f3(
        &_tree_sequence,
        num_sample_sets,
        sample_set_sizes,
        samples_sets.data(),
        num_index_tuples,
        set_indexes,
        num_windows,
        NULL,
        TSK_STAT_SITE,
        &f3
    );
    KASSERT(ret == 0, "Failed to compute the f3.", sfkit::assert::light);

    return f3;
}

double TSKitTreeSequence::f2(SampleSet const& sample_set_1, SampleSet const& sample_set_2) const {
    double                f2;
    constexpr int         num_windows     = 0;
    constexpr int         num_sample_sets = 2;
    std::vector<tsk_id_t> samples_sets;

    for (auto const& sample: sample_set_1) {
        samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
    }
    for (auto const& sample: sample_set_2) {
        samples_sets.push_back(asserting_cast<tsk_id_t>(sample));
    }
    tsk_size_t sample_set_sizes[] = {sample_set_1.popcount(), sample_set_2.popcount()};

    tsk_id_t      set_indexes[]    = {0, 1, 2, 3};
    constexpr int num_index_tuples = 1;

    auto ret = tsk_treeseq_f2(
        &_tree_sequence,
        num_sample_sets,
        sample_set_sizes,
        samples_sets.data(),
        num_index_tuples,
        set_indexes,
        num_windows,
        NULL,
        TSK_STAT_SITE,
        &f2
    );
    KASSERT(ret == 0, "Failed to compute the f2.", sfkit::assert::light);

    return f2;
}

double TSKitTreeSequence::num_segregating_sites() const {
    double                num_seg_sites;
    std::vector<tsk_id_t> samples;
    samples.resize(num_samples());

    std::iota(samples.begin(), samples.end(), 0);
    tsk_size_t sample_set_sizes[] = {num_samples()};

    auto ret = tsk_treeseq_segregating_sites(
        &_tree_sequence,
        1,
        sample_set_sizes,
        samples.data(),
        0,
        NULL,
        TSK_STAT_SITE,
        &num_seg_sites
    );
    KASSERT(ret == 0, "Failed to compute the diversity.", sfkit::assert::light);

    return num_seg_sites;
}

std::vector<double> TSKitTreeSequence::allele_frequency_spectrum() const {
    std::vector<double>   results;
    std::vector<tsk_id_t> samples;

    // tskit also outputs an extra 0 (for unwindowed AFS) in the first element of the results vector.
    results.resize(num_samples() + 1);
    samples.resize(num_samples());
    std::iota(samples.begin(), samples.end(), 0);
    std::array<tsk_size_t, 2> sample_set_sizes = {samples.size(), 0};

    auto ret = tsk_treeseq_allele_frequency_spectrum(
        &_tree_sequence,
        1,
        sample_set_sizes.data(),
        samples.data(),
        0,
        nullptr,
        TSK_STAT_POLARISED,
        // reinterpret_cast<double*>(results.data())
        results.data()
    );
    KASSERT(tskit_noerr(ret), "Failed to compute the allele frequency spectrum", sfkit::assert::light);
    return results;
}

TskMutationView TSKitTreeSequence::mutations() const {
    return std::span(_tree_sequence.site_mutations_mem, num_mutations());
}

std::span<tsk_site_t const> TSKitTreeSequence::sites() const {
    KASSERT(_tree_sequence.tree_sites_mem != nullptr);
    return std::span(_tree_sequence.tree_sites_mem, asserting_cast<size_t>(num_sites()));
}

std::span<double const> const TSKitTreeSequence::breakpoints() const {
    double const* breakpoints_ptr = tsk_treeseq_get_breakpoints(&_tree_sequence);
    return std::span(breakpoints_ptr, num_trees());
}

double TSKitTreeSequence::position_of(tsk_id_t site_id) const {
    tsk_site_t ts_site;
    tsk_treeseq_get_site(&_tree_sequence, site_id, &ts_site);
    return ts_site.position;
}

bool TSKitTreeSequence::tskit_noerr(int ret) const {
    return ret == 0;
}

} // namespace sfkit::tskit
