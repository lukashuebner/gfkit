
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <tdt/graph/common.hpp>
#include <tdt/sequence/GenomicSequence.hpp>

#include "tdt/assertion_levels.hpp"

using namespace Catch::Matchers;

TEST_CASE("Mutation", "[GenomicSequenceStorage]") {
    SiteId const       site_id       = GENERATE(0, 1, 2, 3, 4, 5, 6, 7, 8, 9);
    TreeId const       tree_id       = GENERATE(0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u);
    AllelicState const allelic_state = GENERATE('A', 'C', 'G', 'T');
    AllelicState const parent_state  = GENERATE('A', 'C', 'G', 'T');
    NodeId const       node_id       = GENERATE(1u, 2u, 4u, 5u);

    Mutation mutation(site_id, tree_id, node_id, allelic_state, parent_state);
    CHECK(mutation.site_id() == site_id);
    CHECK(mutation.node_id() == node_id);
    CHECK(mutation.allelic_state() == allelic_state);
}

TEST_CASE("GenomicSequenceStorage Basics", "[GenomicSequenceStorage]") {
    size_t const num_sites_hint    = 10;
    size_t const num_subtrees_hint = 10;

    GenomicSequence sequence(num_sites_hint, num_subtrees_hint);

    // push_back() and get() a few sites
    CHECK(sequence.num_sites() == 0);
    SiteId num_sites = 0;
    for (AllelicState allelic_state: {'A', 'C', 'G', 'T'}) {
        sequence.push_back(allelic_state);
        ++num_sites;
    }
    CHECK(sequence.num_sites() == num_sites);
    for (SiteId site_id = 0; site_id < num_sites; ++site_id) {
        CHECK(sequence.ancestral_state(site_id) == "ACGT"[site_id]);
        CHECK(sequence[site_id] == "ACGT"[site_id]);
    }

    // Check that we can push more sites and subtrees than we hinted.
    for (size_t i = 0; i < num_sites_hint + 1; ++i) {
        sequence.emplace_back('C');
        num_sites++;
        CHECK(sequence[sequence.num_sites() - 1] == 'C');
        CHECK(sequence.num_sites() == num_sites);
    }

    SECTION("Set a few sites vis set()") {
        AllelicState const allelic_state = GENERATE('A', 'C', 'G', 'T');
        for (SiteId site_id = 0; site_id < sequence.num_sites(); ++site_id) {
            sequence.set(site_id, allelic_state);
            CHECK(sequence.ancestral_state(site_id) == allelic_state);
        }
    }

    SECTION("Set a few sites vis operator[]") {
        AllelicState const allelic_state = GENERATE('A', 'C', 'G', 'T');
        for (SiteId site_id = 0; site_id < sequence.num_sites(); ++site_id) {
            sequence[site_id] = allelic_state;
            CHECK(sequence.ancestral_state(site_id) == allelic_state);
            CHECK(sequence[site_id] == allelic_state);
        }
    }
}

TEST_CASE("GenomicSequenceStorage Mutations", "[GenomicSequenceStorage]") {
    GenomicSequence sequence;

    // Adding sites does not affect the number of mutations.
    sequence.push_back('A');
    sequence.push_back('T');
    sequence.push_back('C');
    sequence.push_back('G');
    sequence.emplace_back('G');
    CHECK(sequence.num_mutations() == 0);
    CHECK(sequence.num_sites() == 5);

    // Adding mutations does affect the number of sites.
    SiteId       num_sites    = sequence.num_sites();
    SiteId const site_id      = GENERATE(0, 1, 2, 3, 4, 5, 6, 7, 8, 9);
    TreeId const tree_id      = GENERATE(0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u);
    NodeId const node_id      = GENERATE(1u, 2u, 4u, 5u);
    AllelicState parent_state = 'A';
    for (AllelicState const allelic_state: {'A', 'C', 'G', 'T'}) {
        sequence.emplace_back(site_id, tree_id, node_id, allelic_state, parent_state);
        sequence.push_back(Mutation{site_id, tree_id, node_id, allelic_state, parent_state});
        CHECK(sequence.num_sites() == num_sites);
    }
    CHECK(sequence.num_mutations() == 8);
}
