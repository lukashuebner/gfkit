
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>

#include "sfkit/SuccinctForest.hpp"
#include "sfkit/assertion_levels.hpp"
#include "sfkit/graph/CompressedForest.hpp"
#include "sfkit/load/ForestCompressor.hpp"
#include "sfkit/sequence/AlleleFrequencySpectrum.hpp"
#include "sfkit/tskit.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch;
using namespace ::Catch::Matchers;

// This test case is taken from the tskit test suite (the only test case for the AFS in there that checks values).
TEST_CASE("Tajima's D tskit example", "[Tajima's D]") {
    struct Testcase {
        std::string ts_file;
        double      tskit_tajimas_d;
    };

    // 2, 3, and are multiallelic
    std::vector<Testcase> const testcases = {
        {"data/test-sarafina.trees", -0.45947037060657214},
        {"data/test-scar.trees", -0.4279734893252207},
        {"data/test-shenzi.trees", -0.5179950183615757},
        {"data/test-banzai.trees", -2.2672645402359284},
        {"data/test-ed.trees", -2.672895036702143},
        {"data/test-simba.trees", -2.0469586607448598},
    };
    auto const& testcase        = GENERATE_REF(from_range(testcases));
    auto const& ts_file         = testcase.ts_file;
    auto const& tskit_tajimas_d = testcase.tskit_tajimas_d;
    std::cout << "Testing " << ts_file << " with expected value " << tskit_tajimas_d << std::endl;

    SuccinctForest sequence_forest(ts_file);
    CHECK(Approx(sequence_forest.tajimas_d()).epsilon(1e-6) == tskit_tajimas_d);
}
