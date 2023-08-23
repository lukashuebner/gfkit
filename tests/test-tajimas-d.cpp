
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <kassert/kassert.hpp>

#include "tdt/assertion_levels.hpp"
#include "tdt/graph/compressed-forest.hpp"
#include "tdt/load/forest-compressor.hpp"
#include "tdt/succinct-forest.hpp"
#include "tdt/sequence/allele-frequency-spectrum.hpp"
#include "tdt/tskit.hpp"
#include "tskit-testlib/testlib.hpp"

using namespace ::Catch;
using namespace ::Catch::Matchers;

// This test case is taken from the tskit test suite (the only test case for the AFS in there that checks values).
TEST_CASE("Tajima's D tskit example", "[Diversity]") {
    struct Testcase {
        std::string ts_file;
        double      tskit_tajimas_d;
    };

    // 2, 3, and are multiallelic
    std::vector<Testcase> const testcases = {
        {"data/allele-frequency-spectrum-simple-example-0.trees", -0.45947037060657214},
        {"data/allele-frequency-spectrum-simple-example-1.trees", -0.4279734893252207},
    };
    auto const& testcase        = GENERATE_REF(from_range(testcases));
    auto const& ts_file         = testcase.ts_file;
    auto const& tskit_tajimas_d = testcase.tskit_tajimas_d;
    std::cout << "Testing " << ts_file << " with expected value " << tskit_tajimas_d << std::endl;

    SuccinctForest sequence_forest(ts_file);
    CHECK(Approx(sequence_forest.tajimas_d()).epsilon(1e-6) == tskit_tajimas_d);
}
