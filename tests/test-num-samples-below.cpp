#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <tskit.h>

#include "mocks/TsToSfMappingExtractor.hpp"
#include "tdt/assertion_levels.hpp"
#include "tdt/graph/CompressedForest.hpp"
#include "tdt/load/CompressedForestIO.hpp"
#include "tdt/samples/NumSamplesBelow.hpp"
#include "tdt/samples/SampleSet.hpp"
#include "tdt/sequence/GenomicSequence.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch::Matchers;

TEST_CASE("NumSamplesBelow", "[NumSamplesBelow]") {
    EdgeListGraph dag;
    dag.num_nodes(0);

    SECTION("Empty graph") {
        auto num_samples_below = NumSamplesBelowFactory::build(dag, SampleSet{0});
        CHECK(num_samples_below.num_nodes_in_dag() == 0);
        CHECK(num_samples_below.num_samples_in_dag() == 0);
        CHECK(num_samples_below.num_samples_in_sample_set() == 0);
    }

    SECTION("Graph with a single edge") {
        dag.insert_edge(1, 0);
        dag.insert_leaf(0);
        dag.compute_num_nodes();

        SampleSet empty_samples(dag.num_leaves());
        REQUIRE(empty_samples.overall_num_samples() == 1);
        REQUIRE(empty_samples.popcount() == 0);
        auto empty_num_samples_below = NumSamplesBelowFactory::build(dag, empty_samples);
        CHECK(empty_num_samples_below.num_nodes_in_dag() == 2);
        CHECK(empty_num_samples_below.num_samples_in_dag() == 1);
        CHECK(empty_num_samples_below.num_samples_in_sample_set() == 0);
        CHECK(empty_num_samples_below.num_samples_below(0) == 0);
        CHECK(empty_num_samples_below.num_samples_below(1) == 0);

        SampleSet samples(dag.num_leaves());
        samples.add(0);
        auto num_samples_below = NumSamplesBelowFactory::build(dag, samples);
        CHECK(num_samples_below.num_nodes_in_dag() == 2);
        CHECK(num_samples_below.num_samples_in_dag() == 1);
        CHECK(num_samples_below.num_samples_in_sample_set() == 1);
        CHECK(num_samples_below.num_samples_below(0) == 1);
        CHECK(num_samples_below.num_samples_below(1) == 1);
    }

    SECTION("Graph with three leaves") {
        dag.insert_leaf(0);
        dag.insert_leaf(1);
        dag.insert_leaf(2);
        dag.insert_edge(3, 0);
        dag.insert_edge(3, 1);
        dag.insert_edge(4, 3);
        dag.insert_edge(4, 2);
        dag.compute_num_nodes();

        SECTION("Empty sample set") {
            auto num_samples_below = NumSamplesBelowFactory::build(dag, SampleSet{5});
            CHECK(num_samples_below.num_nodes_in_dag() == 5);
            CHECK(num_samples_below.num_samples_in_dag() == 3);
            CHECK(num_samples_below.num_samples_in_sample_set() == 0);
            CHECK(num_samples_below.num_samples_below(4) == 0);
            CHECK(num_samples_below.num_samples_below(3) == 0);
            CHECK(num_samples_below.num_samples_below(0) == 0);
            CHECK(num_samples_below.num_samples_below(1) == 0);
            CHECK(num_samples_below.num_samples_below(2) == 0);
        }

        SECTION("All samples in sample set") {
            SampleSet samples{5};
            samples.add(0);
            samples.add(1);
            samples.add(2);
            auto num_samples_below = NumSamplesBelowFactory::build(dag, samples);
            CHECK(num_samples_below.num_nodes_in_dag() == 5);
            CHECK(num_samples_below.num_samples_in_dag() == 3);
            CHECK(num_samples_below.num_samples_in_sample_set() == 3);
            CHECK(num_samples_below.num_samples_below(4) == 3);
            CHECK(num_samples_below.num_samples_below(3) == 2);
            CHECK(num_samples_below.num_samples_below(0) == 1);
            CHECK(num_samples_below.num_samples_below(1) == 1);
            CHECK(num_samples_below.num_samples_below(2) == 1);
        }

        SECTION("Some but not all samples in sample set") {
            SampleSet samples{5};
            samples.add(0);
            samples.add(2);
            auto num_samples_below = NumSamplesBelowFactory::build(dag, samples);
            CHECK(num_samples_below.num_nodes_in_dag() == 5);
            CHECK(num_samples_below.num_samples_in_dag() == 3);
            CHECK(num_samples_below.num_samples_in_sample_set() == 2);
            CHECK(num_samples_below.num_samples_below(4) == 2);
            CHECK(num_samples_below.num_samples_below(3) == 1);
            CHECK(num_samples_below.num_samples_below(0) == 1);
            CHECK(num_samples_below.num_samples_below(1) == 0);
            CHECK(num_samples_below.num_samples_below(2) == 1);
        }
    }
}

TEST_CASE("NumSamplesBelow Medium-Sized Example", "[NumSamplesBelow]") {
    //                                      38
    //                                ┏━━━━━━┻━━━━━━━┓
    //                               37              ┃
    //                         ┏━━━━━━┻━━━━━━┓       ┃
    //                        36             ┃       ┃
    //                 ┏━━━━━━━┻━━━━━━━━┓    ┃       ┃
    //                 ┃                ┃   35       ┃
    //                 ┃                ┃   ┏┻┓      ┃
    //                 ┃                ┃   ┃ ┃     34
    //                 ┃                ┃   ┃ ┃    ┏━┻━┓
    //                33                ┃   ┃ ┃    ┃   ┃
    //        ┏━━━━━━━━┻━━━━━━━━┓       ┃   ┃ ┃    ┃   ┃
    //       32                 ┃       ┃   ┃ ┃    ┃   ┃
    //    ┏━━━┻━━━┓             ┃       ┃   ┃ ┃    ┃   ┃
    //   31       ┃             ┃       ┃   ┃ ┃    ┃   ┃
    //  ┏━┻━┓     ┃             ┃       ┃   ┃ ┃    ┃   ┃
    //  ┃   ┃     ┃             ┃       ┃   ┃ ┃   30   ┃
    //  ┃   ┃     ┃             ┃       ┃   ┃ ┃  ┏━┻┓  ┃
    //  ┃   ┃    29             ┃       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┏━━━┻━━━┓         ┃       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃      28         ┃       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┏━━┻━━┓      ┃       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┃     ┃     27       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┃     ┃   ┏━━┻━━┓    ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┃     ┃  26     ┃    ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┃     ┃ ┏━┻━┓   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃   25     ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃  ┏━┻━━┓  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // 24   ┃ ┃  ┃    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // ┏┻┓  ┃ ┃  ┃    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ 23    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┏┻━┓  ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃   ┃   ┃   22   ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃   ┃   ┃  ┏━┻━┓ ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃  21   ┃  ┃   ┃ ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┏━┻┓  ┃  ┃   ┃ ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┃  ┃  ┃ 20   ┃ ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┃  ┃  ┃ ┏┻┓  ┃ ┃ ┃  ┃  ┃  ┃
    // 0 6 10 8 9 11 12 13 4 7 14 19 1 2 18 3 5 15 16 17

    std::vector<SampleId> const leaves{0, 6, 10, 8, 9, 11, 12, 13, 4, 7, 14, 19, 1, 2, 18, 3, 5, 15, 16, 17};
    SampleId const              root = 38;

    EdgeListGraph dag;
    dag.insert_edge(24, 0);
    dag.insert_edge(24, 6);
    dag.insert_edge(31, 24);
    dag.insert_edge(31, 10);
    dag.insert_edge(29, 8);
    dag.insert_edge(23, 9);
    dag.insert_edge(23, 11);
    dag.insert_edge(25, 12);
    dag.insert_edge(25, 23);
    dag.insert_edge(28, 13);
    dag.insert_edge(28, 25);
    dag.insert_edge(29, 28);
    dag.insert_edge(26, 4);
    dag.insert_edge(21, 7);
    dag.insert_edge(21, 14);
    dag.insert_edge(26, 21);
    dag.insert_edge(27, 19);
    dag.insert_edge(27, 26);
    dag.insert_edge(32, 31);
    dag.insert_edge(32, 29);
    dag.insert_edge(33, 32);
    dag.insert_edge(33, 27);
    dag.insert_edge(20, 1);
    dag.insert_edge(20, 2);
    dag.insert_edge(22, 18);
    dag.insert_edge(22, 20);
    dag.insert_edge(36, 33);
    dag.insert_edge(36, 22);
    dag.insert_edge(35, 3);
    dag.insert_edge(35, 5);
    dag.insert_edge(37, 36);
    dag.insert_edge(37, 35);
    dag.insert_edge(30, 15);
    dag.insert_edge(30, 16);
    dag.insert_edge(34, 30);
    dag.insert_edge(34, 17);
    dag.insert_edge(38, 37);
    dag.insert_edge(38, 34);

    for (auto const& leaf: leaves) {
        dag.insert_leaf(leaf);
    }
    dag.insert_root(root);
    dag.compute_num_nodes();

    SECTION("All samples in sample set") {
        SampleSet samples{39};
        samples.add(leaves);

        auto num_samples_below = NumSamplesBelowFactory::build(dag, samples);
        CHECK(num_samples_below.num_nodes_in_dag() == 39);
        CHECK(num_samples_below.num_samples_in_dag() == leaves.size());
        CHECK(num_samples_below.num_samples_in_sample_set() == leaves.size());
        CHECK(num_samples_below.num_samples_below(38) == leaves.size());
        CHECK(num_samples_below.num_samples_below(37) == 17);
        CHECK(num_samples_below.num_samples_below(34) == 3);
        CHECK(num_samples_below.num_samples_below(22) == 3);
        CHECK(num_samples_below.num_samples_below(28) == 4);
        CHECK(num_samples_below.num_samples_below(0) == 1);
        CHECK(num_samples_below.num_samples_below(4) == 1);
    }

    SECTION("Subset of samples in sample set") {
        SampleSet samples{39};
        samples.add(std::vector<SampleId>{6, 10, 8, 9, 11, 12, 13, 4, 7, 14, 19, 3, 5, 15, 16, 17});
        REQUIRE(samples.popcount() == 16);

        auto num_samples_below = NumSamplesBelowFactory::build(dag, samples);
        CHECK(num_samples_below.num_nodes_in_dag() == 39);
        CHECK(num_samples_below.num_samples_in_dag() == leaves.size());
        CHECK(num_samples_below.num_samples_in_sample_set() == 16);
        CHECK(num_samples_below.num_samples_below(38) == samples.popcount());
        CHECK(num_samples_below.num_samples_below(37) == 13);
        CHECK(num_samples_below.num_samples_below(34) == 3);
        CHECK(num_samples_below.num_samples_below(22) == 0);
        CHECK(num_samples_below.num_samples_below(28) == 4);
        CHECK(num_samples_below.num_samples_below(0) == 0);
        CHECK(num_samples_below.num_samples_below(4) == 1);
    }
}

TEST_CASE("NumSamplesBelow Simultaneous Computation of Two Sample Sets Example", "[NumSamplesBelow]") {
    //                                      38
    //                                ┏━━━━━━┻━━━━━━━┓
    //                               37              ┃
    //                         ┏━━━━━━┻━━━━━━┓       ┃
    //                        36             ┃       ┃
    //                 ┏━━━━━━━┻━━━━━━━━┓    ┃       ┃
    //                 ┃                ┃   35       ┃
    //                 ┃                ┃   ┏┻┓      ┃
    //                 ┃                ┃   ┃ ┃     34
    //                 ┃                ┃   ┃ ┃    ┏━┻━┓
    //                33                ┃   ┃ ┃    ┃   ┃
    //        ┏━━━━━━━━┻━━━━━━━━┓       ┃   ┃ ┃    ┃   ┃
    //       32                 ┃       ┃   ┃ ┃    ┃   ┃
    //    ┏━━━┻━━━┓             ┃       ┃   ┃ ┃    ┃   ┃
    //   31       ┃             ┃       ┃   ┃ ┃    ┃   ┃
    //  ┏━┻━┓     ┃             ┃       ┃   ┃ ┃    ┃   ┃
    //  ┃   ┃     ┃             ┃       ┃   ┃ ┃   30   ┃
    //  ┃   ┃     ┃             ┃       ┃   ┃ ┃  ┏━┻┓  ┃
    //  ┃   ┃    29             ┃       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┏━━━┻━━━┓         ┃       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃      28         ┃       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┏━━┻━━┓      ┃       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┃     ┃     27       ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┃     ┃   ┏━━┻━━┓    ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┃     ┃  26     ┃    ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃    ┃     ┃ ┏━┻━┓   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃   25     ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    //  ┃   ┃ ┃  ┏━┻━━┓  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // 24   ┃ ┃  ┃    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // ┏┻┓  ┃ ┃  ┃    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ 23    ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┏┻━┓  ┃  ┃ ┃   ┃   ┃    ┃   ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃   ┃   ┃   22   ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃   ┃   ┃  ┏━┻━┓ ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃  21   ┃  ┃   ┃ ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┏━┻┓  ┃  ┃   ┃ ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┃  ┃  ┃ 20   ┃ ┃ ┃  ┃  ┃  ┃
    // ┃ ┃  ┃ ┃ ┃  ┃  ┃  ┃ ┃ ┃  ┃  ┃ ┏┻┓  ┃ ┃ ┃  ┃  ┃  ┃
    // 0 6 10 8 9 11 12 13 4 7 14 19 1 2 18 3 5 15 16 17

    std::vector<SampleId> const leaves{0, 6, 10, 8, 9, 11, 12, 13, 4, 7, 14, 19, 1, 2, 18, 3, 5, 15, 16, 17};
    SampleId const              root = 38;

    EdgeListGraph dag;
    dag.insert_edge(24, 0);
    dag.insert_edge(24, 6);
    dag.insert_edge(31, 24);
    dag.insert_edge(31, 10);
    dag.insert_edge(29, 8);
    dag.insert_edge(23, 9);
    dag.insert_edge(23, 11);
    dag.insert_edge(25, 12);
    dag.insert_edge(25, 23);
    dag.insert_edge(28, 13);
    dag.insert_edge(28, 25);
    dag.insert_edge(29, 28);
    dag.insert_edge(26, 4);
    dag.insert_edge(21, 7);
    dag.insert_edge(21, 14);
    dag.insert_edge(26, 21);
    dag.insert_edge(27, 19);
    dag.insert_edge(27, 26);
    dag.insert_edge(32, 31);
    dag.insert_edge(32, 29);
    dag.insert_edge(33, 32);
    dag.insert_edge(33, 27);
    dag.insert_edge(20, 1);
    dag.insert_edge(20, 2);
    dag.insert_edge(22, 18);
    dag.insert_edge(22, 20);
    dag.insert_edge(36, 33);
    dag.insert_edge(36, 22);
    dag.insert_edge(35, 3);
    dag.insert_edge(35, 5);
    dag.insert_edge(37, 36);
    dag.insert_edge(37, 35);
    dag.insert_edge(30, 15);
    dag.insert_edge(30, 16);
    dag.insert_edge(34, 30);
    dag.insert_edge(34, 17);
    dag.insert_edge(38, 37);
    dag.insert_edge(38, 34);

    for (auto const& leaf: leaves) {
        dag.insert_leaf(leaf);
    }
    dag.insert_root(root);
    dag.compute_num_nodes();

    SampleSet samples_0{39};
    SampleSet samples_1{39};
    samples_0.add(std::vector<SampleId>{6, 10, 8, 9, 11, 12, 13, 4, 7, 14, 19, 3, 5, 15, 16, 17});
    samples_1.add(std::vector<SampleId>{0, 1, 2, 18});
    REQUIRE(samples_0.popcount() == 16);
    REQUIRE(samples_1.popcount() == 4);

    auto num_samples_below_0 = NumSamplesBelowFactory::build(dag, samples_0);
    auto num_samples_below_1 = NumSamplesBelowFactory::build(dag, samples_1);
    auto [num_samples_below_combined_0, num_samples_below_combined_1] =
        NumSamplesBelowFactory::build(dag, samples_0, samples_1);
    CHECK(num_samples_below_0.num_nodes_in_dag() == num_samples_below_1.num_nodes_in_dag());
    CHECK(num_samples_below_0.num_nodes_in_dag() == num_samples_below_combined_0.num_nodes_in_dag());
    CHECK(num_samples_below_1.num_nodes_in_dag() == num_samples_below_combined_1.num_nodes_in_dag());
    CHECK(num_samples_below_1.num_samples_in_dag() == leaves.size());
    CHECK(num_samples_below_0.num_samples_in_sample_set() == 16);
    CHECK(num_samples_below_1.num_samples_in_sample_set() == 4);
    CHECK(num_samples_below_0.num_samples_below(38) == samples_0.popcount());
    CHECK(num_samples_below_1.num_samples_below(38) == samples_1.popcount());

    for (NodeId node = 0; node < dag.num_nodes(); ++node) {
        CHECK(num_samples_below_0.num_samples_below(node) == num_samples_below_combined_0(node));
        CHECK(num_samples_below_1.num_samples_below(node) == num_samples_below_combined_1(node));
        CHECK(num_samples_below_0.num_samples_below(node) == num_samples_below_0.num_samples_below(node));
        CHECK(num_samples_below_1.num_samples_below(node) == num_samples_below_1.num_samples_below(node));
    }
}

TEST_CASE("NumSamplesBelow Simultaneous Computation of Multiple Sample Sets Simulated", "[NumSamplesBelow]") {
    std::vector<std::string> const ts_files = {
        "data/test-sarafina.trees",
        "data/test-scar.trees",
        "data/test-shenzi.trees",
        "data/test-banzai.trees",
        "data/test-ed.trees",
        "data/test-simba.trees",
    };
    auto const& ts_file = GENERATE_REF(from_range(ts_files));

    TSKitTreeSequence tree_sequence(ts_file);

    ForestCompressor       forest_compressor(tree_sequence);
    GenomicSequenceFactory sequence_factory(tree_sequence);
    CompressedForest       forest = forest_compressor.compress(sequence_factory);

    { // Two sample sets at once
        SampleSet sample_set_0(forest.num_samples());
        SampleSet sample_set_1(forest.num_samples());
        bool      flip = false;
        for (SampleId sample: forest.leaves()) {
            if (flip) {
                sample_set_0.add(sample);
            } else {
                sample_set_1.add(sample);
            }
            flip = !flip;
        }
        REQUIRE(forest.leaves().size() == forest.num_samples());
        REQUIRE(sample_set_0.popcount() + sample_set_1.popcount() == forest.num_samples());

        auto num_samples_below_0_ref = NumSamplesBelowFactory::build(forest.postorder_edges(), sample_set_0);
        auto num_samples_below_1_ref = NumSamplesBelowFactory::build(forest.postorder_edges(), sample_set_1);

        auto [num_samples_below_0, num_samples_below_1] =
            NumSamplesBelowFactory::build(forest.postorder_edges(), sample_set_0, sample_set_1);

        for (NodeId node = 0; node < forest.num_nodes(); ++node) {
            CHECK(num_samples_below_0_ref(node) == num_samples_below_0(node));
            CHECK(num_samples_below_1_ref(node) == num_samples_below_1(node));
            CHECK(num_samples_below_0(node) + num_samples_below_1(node) <= forest.num_samples());
        }

        for (auto& root_node: forest.roots()) {
            CHECK(num_samples_below_0(root_node) + num_samples_below_1(root_node) == forest.num_samples());
        }
    }

    { // Four sample sets at once
        constexpr size_t       num_sample_sets = 4;
        std::vector<SampleSet> sample_sets(num_sample_sets, forest.num_samples());

        size_t idx = 0;
        for (SampleId sample: forest.leaves()) {
            KASSERT(sample < forest.num_samples(), "Sample id out of range", tdt::assert::light);
            sample_sets[idx].add(sample);
            idx = (idx + 1ul) % num_sample_sets;
        }
        REQUIRE(forest.leaves().size() == forest.num_samples());
        REQUIRE(
            sample_sets[0].popcount() + sample_sets[1].popcount() + sample_sets[2].popcount()
                + sample_sets[3].popcount()
            == forest.num_samples()
        );

        auto num_samples_below_0_ref = NumSamplesBelowFactory::build(forest.postorder_edges(), sample_sets[0]);
        auto num_samples_below_1_ref = NumSamplesBelowFactory::build(forest.postorder_edges(), sample_sets[1]);
        auto num_samples_below_2_ref = NumSamplesBelowFactory::build(forest.postorder_edges(), sample_sets[2]);
        auto num_samples_below_3_ref = NumSamplesBelowFactory::build(forest.postorder_edges(), sample_sets[3]);

        auto [num_samples_below_0, num_samples_below_1, num_samples_below_2, num_samples_below_3] =
            NumSamplesBelowFactory::build(
                forest.postorder_edges(),
                sample_sets[0],
                sample_sets[1],
                sample_sets[2],
                sample_sets[3]
            );

        for (NodeId node = 0; node < forest.num_nodes(); ++node) {
            CHECK(num_samples_below_0_ref(node) == num_samples_below_0(node));
            CHECK(num_samples_below_1_ref(node) == num_samples_below_1(node));
            CHECK(num_samples_below_2_ref(node) == num_samples_below_2(node));
            CHECK(num_samples_below_3_ref(node) == num_samples_below_3(node));
            CHECK(
                num_samples_below_0(node) + num_samples_below_1(node) + num_samples_below_2(node)
                    + num_samples_below_3(node)
                <= forest.num_samples()
            );
        }

        for (auto& root_node: forest.roots()) {
            CHECK(
                num_samples_below_0(root_node) + num_samples_below_1(root_node) + num_samples_below_2(root_node)
                    + num_samples_below_3(root_node)
                == forest.num_samples()
            );
        }
    }
}
