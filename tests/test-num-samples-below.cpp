#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>
#include <tskit.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/samples/num-samples-below.hpp"
#include "tdt/samples/sample-set.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch::Matchers;

TEST_CASE("SampleSet", "[NumSamplesBelow]") {
    SampleSet samples{10};

    CHECK(samples.num_nodes_in_dag() == 10);
    CHECK(samples.popcount() == 0);
    // CHECK(samples.build_inverse().overall_num_samples() == samples.overall_num_samples());
    // CHECK(samples.build_inverse().popcount() == 10);

    samples.add(0);
    samples.add(9);
    samples.add(5);

    // SampleSet inverse = samples.build_inverse();
    CHECK(samples.num_nodes_in_dag() == 10);
    CHECK(samples.popcount() == 3);
    // CHECK(inverse.overall_num_samples() == samples.overall_num_samples());
    // CHECK(inverse.popcount() == 7);

    // Check the operator[]
    CHECK(samples[0] == true);
    CHECK(samples[1] == false);
    CHECK(samples[2] == false);
    CHECK(samples[3] == false);
    CHECK(samples[4] == false);
    CHECK(samples[5] == true);
    CHECK(samples[6] == false);
    CHECK(samples[7] == false);
    CHECK(samples[8] == false);
    CHECK(samples[9] == true);

    // CHECK(inverse[0] == false);
    // CHECK(inverse[1] == true);
    // CHECK(inverse[2] == true);
    // CHECK(inverse[3] == true);
    // CHECK(inverse[4] == true);
    // CHECK(inverse[5] == false);
    // CHECK(inverse[6] == true);
    // CHECK(inverse[7] == true);
    // CHECK(inverse[8] == true);
    // CHECK(inverse[9] == false);

    // Check the iterator
    CHECK_THAT(samples, RangeEquals(std::vector<SampleId>{0, 5, 9}));
    // CHECK_THAT(inverse, RangeEquals(std::vector<SampleId>{1, 2, 3, 4, 6, 7, 8}));
}

TEST_CASE("NumSamplesBelow", "[NumSamplesBelow]") {
    EdgeListGraph dag;

    SECTION("Empty graph") {
        NumSamplesBelow num_samples_below(dag, SampleSet{0});
        CHECK(num_samples_below.num_nodes_in_dag() == 0);
        CHECK(num_samples_below.num_samples_in_dag() == 0);
        CHECK(num_samples_below.num_samples_in_sample_set() == 0);
    }

    SECTION("Graph with a single edge") {
        dag.add_edge(0, 1);
        dag.add_leaf(1);
        dag.compute_nodes();

        SampleSet empty_samples(dag.num_leaves());
        REQUIRE(empty_samples.num_nodes_in_dag() == 1);
        REQUIRE(empty_samples.popcount() == 0);
        NumSamplesBelow empty_num_samples_below(dag, empty_samples);
        CHECK(empty_num_samples_below.num_nodes_in_dag() == 2);
        CHECK(empty_num_samples_below.num_samples_in_dag() == 1);
        CHECK(empty_num_samples_below.num_samples_in_sample_set() == 0);
        CHECK(empty_num_samples_below.num_samples_below(0) == 0);
        CHECK(empty_num_samples_below.num_samples_below(1) == 0);

        SampleSet samples(dag.num_nodes());
        samples.add(1);
        NumSamplesBelow num_samples_below(dag, samples);
        CHECK(num_samples_below.num_nodes_in_dag() == 2);
        CHECK(num_samples_below.num_samples_in_dag() == 1);
        CHECK(num_samples_below.num_samples_in_sample_set() == 1);
        CHECK(num_samples_below.num_samples_below(0) == 1);
        CHECK(num_samples_below.num_samples_below(1) == 1);
    }

    SECTION("Graph with three leaves") {
        dag.add_leaf(2);
        dag.add_leaf(3);
        dag.add_leaf(4);
        dag.add_edge(1, 2);
        dag.add_edge(1, 3);
        dag.add_edge(0, 4);
        dag.add_edge(0, 1);
        dag.compute_nodes();

        SECTION("Empty sample set") {
            NumSamplesBelow num_samples_below(dag, SampleSet{5});
            CHECK(num_samples_below.num_nodes_in_dag() == 5);
            CHECK(num_samples_below.num_samples_in_dag() == 3);
            CHECK(num_samples_below.num_samples_in_sample_set() == 0);
            CHECK(num_samples_below.num_samples_below(0) == 0);
            CHECK(num_samples_below.num_samples_below(1) == 0);
            CHECK(num_samples_below.num_samples_below(2) == 0);
            CHECK(num_samples_below.num_samples_below(3) == 0);
            CHECK(num_samples_below.num_samples_below(4) == 0);
        }

        SECTION("All samples in sample set") {
            SampleSet samples{5};
            samples.add(2);
            samples.add(3);
            samples.add(4);
            NumSamplesBelow num_samples_below(dag, samples);
            CHECK(num_samples_below.num_nodes_in_dag() == 5);
            CHECK(num_samples_below.num_samples_in_dag() == 3);
            CHECK(num_samples_below.num_samples_in_sample_set() == 3);
            CHECK(num_samples_below.num_samples_below(0) == 3);
            CHECK(num_samples_below.num_samples_below(1) == 2);
            CHECK(num_samples_below.num_samples_below(2) == 1);
            CHECK(num_samples_below.num_samples_below(3) == 1);
            CHECK(num_samples_below.num_samples_below(4) == 1);
        }

        SECTION("Some but not all samples in sample set") {
            SampleSet samples{5};
            samples.add(2);
            samples.add(4);
            NumSamplesBelow num_samples_below(dag, samples);
            CHECK(num_samples_below.num_nodes_in_dag() == 5);
            CHECK(num_samples_below.num_samples_in_dag() == 3);
            CHECK(num_samples_below.num_samples_in_sample_set() == 2);
            CHECK(num_samples_below.num_samples_below(0) == 2);
            CHECK(num_samples_below.num_samples_below(1) == 1);
            CHECK(num_samples_below.num_samples_below(2) == 1);
            CHECK(num_samples_below.num_samples_below(3) == 0);
            CHECK(num_samples_below.num_samples_below(4) == 1);
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
    dag.add_edge(24, 0);
    dag.add_edge(24, 6);
    dag.add_edge(31, 24);
    dag.add_edge(31, 10);
    dag.add_edge(29, 8);
    dag.add_edge(23, 9);
    dag.add_edge(23, 11);
    dag.add_edge(25, 12);
    dag.add_edge(25, 23);
    dag.add_edge(28, 13);
    dag.add_edge(28, 25);
    dag.add_edge(29, 28);
    dag.add_edge(26, 4);
    dag.add_edge(21, 7);
    dag.add_edge(21, 14);
    dag.add_edge(26, 21);
    dag.add_edge(27, 19);
    dag.add_edge(27, 26);
    dag.add_edge(32, 31);
    dag.add_edge(32, 29);
    dag.add_edge(33, 32);
    dag.add_edge(33, 27);
    dag.add_edge(20, 1);
    dag.add_edge(20, 2);
    dag.add_edge(22, 18);
    dag.add_edge(22, 20);
    dag.add_edge(36, 33);
    dag.add_edge(36, 22);
    dag.add_edge(35, 3);
    dag.add_edge(35, 5);
    dag.add_edge(37, 36);
    dag.add_edge(37, 35);
    dag.add_edge(30, 15);
    dag.add_edge(30, 16);
    dag.add_edge(34, 30);
    dag.add_edge(34, 17);
    dag.add_edge(38, 37);
    dag.add_edge(38, 34);

    for (auto const& leaf: leaves) {
        dag.add_leaf(leaf);
    }
    dag.add_root(root);
    dag.compute_nodes();

    SECTION("All samples in sample set") {
        SampleSet samples{39};
        samples.add(leaves);

        NumSamplesBelow num_samples_below(dag, samples);
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

        NumSamplesBelow num_samples_below(dag, samples);
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
