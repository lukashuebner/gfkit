#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <tskit.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/graph/compressed-forest.hpp"
#include "tdt/load/forest-compressor.hpp"
#include "tdt/sequence-forest.hpp"
#include "tdt/tskit.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch::Matchers;

// This test case is taken from the tskit test suite (the only test case for the AFS in there that checks values).
TEST_CASE("AlleleFrequencies example multi tree no back no recurrent", "[AlleleFrequencies]") {
    tsk_treeseq_t tskit_tree_sequence;

    tsk_treeseq_from_text(
        &tskit_tree_sequence,
        10,
        paper_ex_nodes,
        paper_ex_edges,
        NULL,
        paper_ex_sites,
        paper_ex_mutations,
        paper_ex_individuals,
        NULL,
        0
    );

    TSKitTreeSequence tdt_tree_sequence(tskit_tree_sequence); // Takes ownership
    ForestCompressor  forest_compressor(tdt_tree_sequence);
    CompressedForest  compressed_forest = forest_compressor.compress();

    GenomicSequenceStorage sequence_store(tdt_tree_sequence, forest_compressor);
    SequenceForest         sequence_forest(compressed_forest, sequence_store);

    // .allele_frequencies() returns the number of samples in the ancestral state.
    std::vector<uint64_t> expected = {3, 3, 1};
    CHECK_THAT(sequence_forest.allele_frequencies(), RangeEquals(expected));

    // Do not free the tskit tree sequence, as we transferred ownershop to  tdt_tree_sequence now.
    // tsk_treeseq_free(&tskit_tree_sequence);
}
