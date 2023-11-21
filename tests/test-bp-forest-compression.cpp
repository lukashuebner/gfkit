#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_quantifiers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>

#include "mocks/TsToSfMappingExtractor.hpp"
#include "sfkit/SuccinctForest.hpp"
#include "sfkit/assertion_levels.hpp"
#include "sfkit/graph/AdjacencyArrayGraph.hpp"
#include "sfkit/graph/EdgeListGraph.hpp"
#include "sfkit/graph/common.hpp"
#include "sfkit/load/BPForestCompressor.hpp"
#include "sfkit/load/ForestCompressor.hpp"
#include "sfkit/utils/concepts.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace Catch::Matchers;

TEST_CASE("BP Forest Compression Example I", "[BPForestCompression]") {
    /*          6          */
    /*         / \         */
    /*        /   \        */
    /*       /     \       */
    /*      /       5      */
    /*     4       / \     */
    /*    / \     /   \    */
    /*   0   1   2     3   */
    tsk_treeseq_t tskit_tree_sequence;

    tsk_treeseq_from_text(
        &tskit_tree_sequence,
        1,
        single_tree_multi_derived_states_nodes,
        single_tree_multi_derived_states_edges,
        NULL,
        single_tree_multi_derived_states_sites,
        single_tree_multi_derived_states_mutations,
        NULL,
        NULL,
        0
    );

    TSKitTreeSequence     tree_sequence(tskit_tree_sequence); // Takes ownership
    Ts2SfMappingExtractor ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    auto const            bp_forest = BPForestCompressor(tree_sequence).compress(ts_2_sf_node);
    CHECK(ts_2_sf_node.finalize_called());
    CHECK(ts_2_sf_node.process_mutations_callcnt() == 1);

    Ts2SfMappingExtractor ts_2_sf_node_2(tree_sequence.num_trees(), tree_sequence.num_nodes());
    auto const            non_bp_forest = ForestCompressor(tree_sequence).compress(ts_2_sf_node);

    CHECK(bp_forest.num_trees() == 1);
    CHECK(bp_forest.num_samples() == 4);
    CHECK(bp_forest.num_leaves() == bp_forest.num_samples());
    CHECK(bp_forest.leaves().size() == bp_forest.num_leaves());
    CHECK(bp_forest.num_unique_subtrees() == 7);
    CHECK(bp_forest.num_nodes() == bp_forest.num_unique_subtrees());

    CHECK_THAT(bp_forest.is_reference(), NoneTrue());
    CHECK_THAT(
        bp_forest.balanced_parenthesis(),
        RangeEquals(std::vector<bool>{1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0})
    );
    CHECK_THAT(bp_forest.is_leaf(), RangeEquals(std::vector<bool>{0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0}));
    CHECK_THAT(bp_forest.leaves(), RangeEquals(std::vector<NodeId>{0, 1, 2, 3}));
    CHECK(bp_forest.references().size() == 0);

    SampleSet const all_samples       = bp_forest.all_samples();
    auto            num_samples_below = NumSamplesBelowFactory::build(bp_forest, all_samples);
    auto num_samples_below_reference  = NumSamplesBelowFactory::build(non_bp_forest.postorder_edges(), all_samples);
    // TODO (also below)
    //! This depends on the node_ids being the same in both compression, which they generally aren't. See the zazu test
    //! for code which does not assume this.
    for (NodeId node_id = 0; node_id < bp_forest.num_nodes(); ++node_id) {
        CHECK(num_samples_below(node_id) == num_samples_below_reference(node_id));
    }

    // Do not free the tskit tree sequence, as we transferred the ownershop to sfkit_tree_sequence.
    // tsk_treeseq_free(&tskit_tree_sequence);
}

TEST_CASE("BP Forest Compression Example II", "[BPForestCompression]") {
    /*
    0.25┊     8   ┊         ┊         ┊
        ┊   ┏━┻━┓ ┊         ┊         ┊
    0.20┊   ┃   ┃ ┊         ┊   7     ┊
        ┊   ┃   ┃ ┊         ┊ ┏━┻━┓   ┊
    0.17┊   6   ┃ ┊   6     ┊ ┃   ┃   ┊
        ┊ ┏━┻┓  ┃ ┊ ┏━┻━┓   ┊ ┃   ┃   ┊
    0.09┊ ┃  5  ┃ ┊ ┃   5   ┊ ┃   5   ┊
        ┊ ┃ ┏┻┓ ┃ ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊
    0.07┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊
        ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊
    0.00┊ 0 1 3 2 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊
    */

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
    constexpr TreeId   num_trees           = 3;
    constexpr SampleId num_samples         = 4;
    constexpr NodeId   num_unique_subtrees = 10;
    constexpr NodeId   num_references      = 5;

    TSKitTreeSequence     tree_sequence(tskit_tree_sequence); // Takes ownership
    Ts2SfMappingExtractor ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    auto const            bp_forest = BPForestCompressor(tree_sequence).compress(ts_2_sf_node);
    CHECK(ts_2_sf_node.finalize_called());
    CHECK(ts_2_sf_node.process_mutations_callcnt() == num_trees);

    Ts2SfMappingExtractor ts_2_sf_node_2(tree_sequence.num_trees(), tree_sequence.num_nodes());
    auto const            non_bp_forest = ForestCompressor(tree_sequence).compress(ts_2_sf_node);

    CHECK(bp_forest.num_trees() == num_trees);
    CHECK(bp_forest.num_samples() == num_samples);
    CHECK(bp_forest.num_leaves() == bp_forest.num_samples());
    CHECK(bp_forest.leaves().size() == bp_forest.num_leaves());
    CHECK(bp_forest.num_unique_subtrees() == num_unique_subtrees);
    CHECK(bp_forest.num_nodes() == bp_forest.num_unique_subtrees());

    CHECK_THAT(
        bp_forest.balanced_parenthesis(),
        RangeEquals(std::vector<bool>{
            1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, // tree 0
            1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, // tree 1
            1, 0                                      // tree 2
        })
    );

    CHECK_THAT(
        bp_forest.is_reference(),
        RangeEquals(std::vector<bool>{
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // tree 0
            0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, // tree 1
            1, 1                                      // tree 2
        })
    );

    // CHECK_THAT(
    //     bp_forest.is_leaf(),
    //     RangeEquals(std::vector<bool>{
    //         0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, // tree 0
    //         0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, // tree 1
    //         0, 0                                      // tree 2
    //     })
    // );

    CHECK_THAT(
        bp_forest.is_leaf(),
        RangeEquals(std::vector<bool>{
            0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, // tree 0
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // tree 0
            0, 0                                      // tree 2
        })
    );

    // Leaves in the second tree are references to the leaves in the first tree. The order of the samples is not unique,
    // as we don't specify the order in which the children are processed.
    // CHECK_THAT(bp_forest.leaves(), RangeEquals(std::vector<NodeId>{0, 1, 3, 2, 0, 1, 2, 3}));
    CHECK_THAT(bp_forest.leaves(), RangeEquals(std::vector<NodeId>{2, 0, 1, 3}));

    // Each of the four samples is referenced in the second tree. The whole second tree is referenced in the third tree.
    // Each reference consists of two elements: <start, length>.
    CHECK(bp_forest.references().size() == num_references);

    SampleSet const all_samples       = bp_forest.all_samples();
    auto const      num_samples_below = NumSamplesBelowFactory::build(bp_forest, all_samples);
    auto const      num_samples_below_reference =
        NumSamplesBelowFactory::build(non_bp_forest.postorder_edges(), all_samples);
    //! This depends on the node_ids being the same in both compression, which they generally aren't. See the zazu test
    //! for code which does not assume this.
    for (NodeId node_id = 0; node_id < bp_forest.num_nodes(); ++node_id) {
        CHECK(num_samples_below(node_id) == num_samples_below_reference(node_id));
    }
}

TEST_CASE("BP Forest Compression Zazu", "[BPForestCompression]") {
    // test-zazu.trees
    // 76.16┊         ┊         ┊         ┊         ┊   9     ┊
    //      ┊         ┊         ┊         ┊         ┊ ┏━┻━┓   ┊
    // 50.34┊   8     ┊         ┊   8     ┊    8    ┊ ┃   8   ┊
    //      ┊ ┏━┻━┓   ┊         ┊ ┏━┻━┓   ┊  ┏━┻━┓  ┊ ┃  ┏┻━┓ ┊
    // 27.71┊ ┃   ┃   ┊         ┊ ┃   ┃   ┊  7   ┃  ┊ ┃  ┃  ┃ ┊
    //      ┊ ┃   ┃   ┊         ┊ ┃   ┃   ┊ ┏┻┓  ┃  ┊ ┃  ┃  ┃ ┊
    // 21.39┊ ┃   ┃   ┊   6     ┊ ┃   ┃   ┊ ┃ ┃  ┃  ┊ ┃  ┃  ┃ ┊
    //      ┊ ┃   ┃   ┊ ┏━┻━┓   ┊ ┃   ┃   ┊ ┃ ┃  ┃  ┊ ┃  ┃  ┃ ┊
    // 11.54┊ ┃   5   ┊ ┃   5   ┊ ┃   5   ┊ ┃ ┃  5  ┊ ┃  5  ┃ ┊
    //      ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┏┻┓ ┃ ┊
    // 5.30 ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊
    //      ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┊
    // 0.00 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊ 0 3 1 2 ┊ 0 1 2 3 ┊
    //      0         2        11        13        19        20
    constexpr TreeId   num_trees           = 5;
    constexpr SampleId num_samples         = 4;
    constexpr NodeId   num_unique_subtrees = 12;
    constexpr NodeId   num_references      = 9;

    TSKitTreeSequence     tree_sequence("data/test-zazu.trees");
    Ts2SfMappingExtractor bp_ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    auto const            bp_forest = BPForestCompressor(tree_sequence).compress(bp_ts_2_sf_node);
    CHECK(bp_ts_2_sf_node.finalize_called());
    CHECK(bp_ts_2_sf_node.process_mutations_callcnt() == num_trees);

    Ts2SfMappingExtractor ref_ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    auto const            non_bp_forest = ForestCompressor(tree_sequence).compress(ref_ts_2_sf_node);
    CHECK(ref_ts_2_sf_node.finalize_called());
    CHECK(ref_ts_2_sf_node.process_mutations_callcnt() == bp_ts_2_sf_node.process_mutations_callcnt());

    CHECK(bp_forest.num_trees() == num_trees);
    CHECK(bp_forest.num_samples() == num_samples);
    CHECK(bp_forest.num_leaves() == bp_forest.num_samples());
    CHECK(bp_forest.leaves().size() == bp_forest.num_leaves());
    CHECK(bp_forest.num_unique_subtrees() == num_unique_subtrees);
    CHECK(bp_forest.num_nodes() == bp_forest.num_unique_subtrees());

    CHECK_THAT(
        bp_forest.balanced_parenthesis(),
        RangeEquals(std::vector<bool>{
            1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, // tree 0
            1, 0,                                     // tree 1
            1, 0,                                     // tree 2
            1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, // tree 3
            1, 1, 0, 1, 1, 0, 1, 0, 0, 0              // tree 4
        })
    );

    CHECK_THAT(
        bp_forest.is_reference(),
        RangeEquals(std::vector<bool>{
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // tree 0
            1, 1,                                     // tree 1
            1, 1,                                     // tree 2
            0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, // tree 3
            0, 1, 1, 0, 1, 1, 1, 1, 0, 0              // tree 4
        })
    );

    // CHECK_THAT(
    //     bp_forest.is_leaf(),
    //     RangeEquals(std::vector<bool>{
    //         0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, // tree 0
    //         0, 0,                                     // tree 1
    //         0, 0,                                     // tree 2
    //         0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, // tree 3
    //         0, 1, 1, 0, 0, 0, 1, 1, 0, 0              // tree 4
    //     })
    // );

    CHECK_THAT(
        bp_forest.is_leaf(),
        RangeEquals(std::vector<bool>{
            0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, // tree 0
            0, 0,                                     // tree 1
            0, 0,                                     // tree 2
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // tree 3
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0              // tree 4
        })
    );

    // Leaves in the second tree are references to the leaves in the first tree. The order of the samples is not unique,
    // as we don't specify the order in which the children are processed.
    CHECK_THAT(bp_forest.leaves(), RangeEquals(std::vector<NodeId>{0, 1, 2, 3}));

    // Each of the four samples is referenced in the second tree. The whole second tree is referenced in the third tree.
    // Each reference consists of two elements: <start, length>.
    CHECK(bp_forest.references().size() == num_references);

    SampleSet const all_samples           = bp_forest.all_samples();
    auto const      bp_num_samples_below  = NumSamplesBelowFactory::build(bp_forest, all_samples);
    auto const      ref_num_samples_below = NumSamplesBelowFactory::build(non_bp_forest.postorder_edges(), all_samples);

    TSKitTree tree{tree_sequence};
    TreeId    tree_id = 0;
    for (tree.first(); tree.is_valid(); tree.next()) {
        KASSERT(tree_id < bp_forest.num_trees(), "Tree id out of range", sfkit::assert::light);
        for (tsk_id_t ts_node_id: tree.postorder()) {
            if (!tree.is_root(ts_node_id)) { // The root nodes are not mapped because they're never referenced
                NodeId const ref_node_id = ref_ts_2_sf_node(tree_id, ts_node_id);
                NodeId const bp_node_id  = bp_ts_2_sf_node(tree_id, ts_node_id);
                if (bp_num_samples_below(bp_node_id) != ref_num_samples_below(ref_node_id)) {
                    std::cout << "Tree: " << tree_id << " Node: " << ts_node_id << std::endl;
                    std::cout << "BP: " << bp_num_samples_below(bp_node_id)
                              << " REF: " << ref_num_samples_below(ref_node_id) << std::endl;
                }
                CHECK(bp_num_samples_below(bp_node_id) == ref_num_samples_below(ref_node_id));
            }
        }
        ++tree_id;
    }
}

TEST_CASE("Compare BP-based compression to reference implementation", "[BPForestCompresion]") {
    std::vector<std::string> const ts_files = {
        "data/test-sarafina.trees",
        "data/test-scar.trees",
        "data/test-shenzi.trees",
        "data/test-banzai.trees",
        "data/test-ed.trees",
        "data/test-simba.trees",
        "data/test-zazu.trees",
        "data/test-pumbaa.trees",
    };
    auto const& ts_file = GENERATE_REF(from_range(ts_files));

    TSKitTreeSequence tree_sequence(ts_file);

    ForestCompressor ref_forest_compressor(tree_sequence);
    // GenomicSequenceFactory ref_sequence_factory(tree_sequence);
    Ts2SfMappingExtractor ref_ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    CompressedForest      ref_forest = ref_forest_compressor.compress(ref_ts_2_sf_node);

    BPForestCompressor bp_forest_compressor(tree_sequence);
    //    GenomicSequenceFactory bp_sequence_factory(tree_sequence);
    Ts2SfMappingExtractor bp_ts_2_sf_node(tree_sequence.num_trees(), tree_sequence.num_nodes());
    BPCompressedForest    bp_forest = bp_forest_compressor.compress(bp_ts_2_sf_node);

    CHECK(bp_forest.num_samples() == ref_forest.num_samples());
    CHECK(bp_forest.num_leaves() == ref_forest.num_leaves());
    CHECK(bp_forest.num_trees() == ref_forest.num_trees());
    CHECK(bp_forest.num_nodes() == ref_forest.num_nodes());
    CHECK(bp_forest.num_unique_subtrees() == ref_forest.num_unique_subtrees());

    { // Single sample set
        SampleSet  all_samples{bp_forest.all_samples()};
        auto const bp_num_samples_below  = NumSamplesBelowFactory::build(bp_forest, all_samples);
        auto const ref_num_samples_below = NumSamplesBelowFactory::build(ref_forest.postorder_edges(), all_samples);

        // TODO Abstract to function
        TSKitTree tree(tree_sequence);
        TreeId    tree_id = 0;
        for (tree.first(); tree.is_valid(); tree.next()) {
            KASSERT(tree_id < bp_forest.num_trees(), "Tree id out of range", sfkit::assert::light);
            for (tsk_id_t ts_node_id: tree.postorder()) {
                if (!tree.is_root(ts_node_id)) { // The root nodes are not mapped because they're never referenced
                    NodeId const ref_node_id = ref_ts_2_sf_node(tree_id, ts_node_id);
                    NodeId const bp_node_id  = bp_ts_2_sf_node(tree_id, ts_node_id);
                    if (bp_num_samples_below(bp_node_id) != ref_num_samples_below(ref_node_id)) {
                        // TODO Remove this output
                        std::cout << "Tree: " << tree_id << " Node: " << ts_node_id << std::endl;
                        std::cout << "BP: " << bp_num_samples_below(bp_node_id)
                                  << " REF: " << ref_num_samples_below(ref_node_id) << std::endl
                                  << "----------" << std::endl;
                    }
                    CHECK(bp_num_samples_below(bp_node_id) == ref_num_samples_below(ref_node_id));
                }
            }
            ++tree_id;
        }
    }

    // TODO
    // { // Two sample sets at once
    //     SampleSet sample_set_0(bp_forest.num_samples());
    //     SampleSet sample_set_1(bp_forest.num_samples());
    //     bool      flip = false;
    //     for (SampleId sample: bp_forest.leaves()) {
    //         if (flip) {
    //             sample_set_0.add(sample);
    //         } else {
    //             sample_set_1.add(sample);
    //         }
    //         flip = !flip;
    //     }
    //     REQUIRE(bp_forest.leaves().size() == bp_forest.num_samples());
    //     REQUIRE(sample_set_0.popcount() + sample_set_1.popcount() == bp_forest.num_samples());

    //     auto num_samples_below_0_ref = NumSamplesBelowFactory::build(forest.postorder_edges(), sample_set_0);
    //     auto num_samples_below_1_ref = NumSamplesBelowFactory::build(forest.postorder_edges(), sample_set_1);

    //     auto [num_samples_below_0, num_samples_below_1] =
    //         NumSamplesBelowFactory::build(forest.postorder_edges(), sample_set_0, sample_set_1);

    //     for (NodeId node = 0; node < forest.num_nodes(); ++node) {
    //         CHECK(num_samples_below_0_ref(node) == num_samples_below_0(node));
    //         CHECK(num_samples_below_1_ref(node) == num_samples_below_1(node));
    //         CHECK(num_samples_below_0(node) + num_samples_below_1(node) <= forest.num_samples());
    //     }

    //     for (auto& root_node: forest.roots()) {
    //         CHECK(num_samples_below_0(root_node) + num_samples_below_1(root_node) == forest.num_samples());
    //     }
    // }

    // { // Four sample sets at once
    //     constexpr size_t       num_sample_sets = 4;
    //     std::vector<SampleSet> sample_sets(num_sample_sets, forest.num_samples());

    //     size_t idx = 0;
    //     for (SampleId sample: forest.leaves()) {
    //         KASSERT(sample < forest.num_samples(), "Sample id out of range", sfkit::assert::light);
    //         sample_sets[idx].add(sample);
    //         idx = (idx + 1ul) % num_sample_sets;
    //     }
    //     REQUIRE(forest.leaves().size() == forest.num_samples());
    //     REQUIRE(
    //         sample_sets[0].popcount() + sample_sets[1].popcount() + sample_sets[2].popcount()
    //             + sample_sets[3].popcount()
    //         == forest.num_samples()
    //     );

    //     auto num_samples_below_0_ref = NumSamplesBelowFactory::build(forest.postorder_edges(), sample_sets[0]);
    //     auto num_samples_below_1_ref = NumSamplesBelowFactory::build(forest.postorder_edges(), sample_sets[1]);
    //     auto num_samples_below_2_ref = NumSamplesBelowFactory::build(forest.postorder_edges(), sample_sets[2]);
    //     auto num_samples_below_3_ref = NumSamplesBelowFactory::build(forest.postorder_edges(), sample_sets[3]);

    //     auto [num_samples_below_0, num_samples_below_1, num_samples_below_2, num_samples_below_3] =
    //         NumSamplesBelowFactory::build(
    //             forest.postorder_edges(),
    //             sample_sets[0],
    //             sample_sets[1],
    //             sample_sets[2],
    //             sample_sets[3]
    //         );

    //     for (NodeId node = 0; node < forest.num_nodes(); ++node) {
    //         CHECK(num_samples_below_0_ref(node) == num_samples_below_0(node));
    //         CHECK(num_samples_below_1_ref(node) == num_samples_below_1(node));
    //         CHECK(num_samples_below_2_ref(node) == num_samples_below_2(node));
    //         CHECK(num_samples_below_3_ref(node) == num_samples_below_3(node));
    //         CHECK(
    //             num_samples_below_0(node) + num_samples_below_1(node) + num_samples_below_2(node)
    //                 + num_samples_below_3(node)
    //             <= forest.num_samples()
    //         );
    //     }

    //     for (auto& root_node: forest.roots()) {
    //         CHECK(
    //             num_samples_below_0(root_node) + num_samples_below_1(root_node) + num_samples_below_2(root_node)
    //                 + num_samples_below_3(root_node)
    //             == forest.num_samples()
    //         );
    //     }
    // }
}
