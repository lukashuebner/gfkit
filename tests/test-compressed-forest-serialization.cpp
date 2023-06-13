
#include <filesystem>
#include <fstream>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>

#include "cereal/archives/binary.hpp"
#include "cereal/cereal.hpp"
#include "tdt/assertion_levels.hpp"
#include "tdt/graph/edge-list-graph.hpp"
#include "tdt/load/compressed-forest-serialization.hpp"
#include "tdt/load/forest-compressor.hpp"
#include "tdt/sequence/allele-frequency-spectrum.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"
#include "tdt/tskit.hpp"

using namespace ::Catch::Matchers;

std::string const ARCHIVE_FILE_NAME = "tmp-test-f5f515340fa29c848db5ed746253f571c6c791bb.forest";

TEST_CASE("CompressedForest/GenomicSequenceStorage Serialization", "[Serialization]") {
    std::vector<std::string> const ts_files = {
        "data/allele-frequency-spectrum-simple-example-0.trees",
        "data/allele-frequency-spectrum-simple-example-1.trees",
        "data/allele-frequency-spectrum-simple-example-2.trees",
        "data/allele-frequency-spectrum-simple-example-3.trees",
        "data/allele-frequency-spectrum-simple-example-4.trees",
        "data/allele-frequency-spectrum-simple-example-6.trees",
    };
    auto const& ts_file = GENERATE_REF(from_range(ts_files));

    TSKitTreeSequence tree_sequence(ts_file);

    auto             reference_afs = tree_sequence.allele_frequency_spectrum();
    ForestCompressor forest_compressor(tree_sequence);
    CompressedForest compressed_forest = forest_compressor.compress();
    compressed_forest.compute_num_samples_below();

    GenomicSequenceStorage sequence_store(tree_sequence, forest_compressor);

    // Serialize and deserialize the compressed forest and genome sequence storage
    CompressedForestIO::save(ARCHIVE_FILE_NAME, compressed_forest, sequence_store);

    CompressedForest       compressed_forest_deserialized;
    GenomicSequenceStorage sequence_store_deserialized;
    CompressedForestIO::load(ARCHIVE_FILE_NAME, compressed_forest_deserialized, sequence_store_deserialized);

    std::filesystem::remove(ARCHIVE_FILE_NAME);

    // CompressedForest
    CHECK(compressed_forest_deserialized.num_nodes() == compressed_forest.num_nodes());
    CHECK(compressed_forest_deserialized.num_samples() == compressed_forest.num_samples());
    CHECK(compressed_forest_deserialized.num_trees() == compressed_forest.num_trees());
    CHECK(compressed_forest_deserialized.num_edges() == compressed_forest.num_edges());
    CHECK(compressed_forest_deserialized.num_roots() == compressed_forest.num_roots());
    CHECK(compressed_forest_deserialized.roots() == compressed_forest.roots());

    for (NodeId node_id = 0; node_id < compressed_forest.num_nodes(); ++node_id) {
        CHECK(compressed_forest_deserialized.is_sample(node_id) == compressed_forest.is_sample(node_id));
        CHECK(
            compressed_forest_deserialized.num_samples_below(node_id) == compressed_forest.num_samples_below(node_id)
        );
    }

    auto cf_dag              = compressed_forest_deserialized.postorder_edges();
    auto cf_dag_deserialized = compressed_forest_deserialized.postorder_edges();
    CHECK(cf_dag_deserialized.num_nodes() == cf_dag.num_nodes());
    CHECK(cf_dag_deserialized.num_edges() == cf_dag.num_edges());
    CHECK(cf_dag_deserialized.num_roots() == cf_dag.num_roots());
    CHECK(cf_dag_deserialized.num_trees() == cf_dag.num_trees());
    CHECK(cf_dag_deserialized.roots() == cf_dag.roots());
    CHECK(cf_dag_deserialized.leaves() == cf_dag.leaves());
    CHECK(cf_dag_deserialized.nodes() == cf_dag.nodes());
    CHECK(cf_dag_deserialized.directed() == cf_dag.directed());
    CHECK(cf_dag_deserialized.traversal_order() == cf_dag.traversal_order());
    CHECK(cf_dag_deserialized.is_inorder() == cf_dag.is_inorder());
    CHECK(cf_dag_deserialized.is_postorder() == cf_dag.is_postorder());
    CHECK(cf_dag_deserialized.is_preorder() == cf_dag.is_preorder());
    CHECK(cf_dag_deserialized.is_levelorder() == cf_dag.is_levelorder());
    CHECK(cf_dag_deserialized.is_unordered() == cf_dag.is_unordered());
    CHECK(cf_dag_deserialized.check_postorder() == cf_dag.check_postorder());
    CHECK_THAT(cf_dag_deserialized, RangeEquals(cf_dag));

    // GenomicSequenceStore
    CHECK(sequence_store_deserialized.num_sites() == sequence_store.num_sites());
    CHECK(sequence_store_deserialized.num_mutations() == sequence_store.num_mutations());

    for (SiteId site_id = 0; site_id < sequence_store.num_sites(); ++site_id) {
        CHECK(sequence_store_deserialized.ancestral_state(site_id) == sequence_store.ancestral_state(site_id));
        CHECK(sequence_store_deserialized[site_id] == sequence_store[site_id]);
        CHECK_THAT(
            sequence_store_deserialized.mutations_at_site(site_id),
            RangeEquals(sequence_store.mutations_at_site(site_id))
        );
    }

    for (size_t mutation_id = 0; mutation_id < sequence_store.num_mutations(); ++mutation_id) {
        CHECK(sequence_store_deserialized.mutation_by_id(mutation_id) == sequence_store.mutation_by_id(mutation_id));
    }

    CHECK(sequence_store_deserialized.mutations_are_sorted_by_site() == sequence_store.mutations_are_sorted_by_site());

    // High-level operations
    AlleleFrequencySpectrum<PerfectDNAHasher> afs(sequence_store, compressed_forest);
    AlleleFrequencySpectrum<PerfectDNAHasher> afs_from_deserialized(
        sequence_store_deserialized,
        compressed_forest_deserialized
    );

    CHECK(afs.num_samples() == afs_from_deserialized.num_samples());
    CHECK(afs.num_sites() == afs_from_deserialized.num_sites());
    CHECK_THAT(afs, RangeEquals(afs_from_deserialized));

    // tskit stores 0's in the lowest and highest position of the AFS for unwindowed AFS, we don't.
    CHECK_THAT(
        std::span(afs).subspan(1, afs.num_samples() - 1),
        RangeEquals(std::span(reference_afs).subspan(1, afs.num_samples() - 1))
    );
}
