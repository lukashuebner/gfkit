list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(Catch)
include(EasyCatch)

# Registering tests:
register_test(test-edge-list-graph FILES test-edge-list-graph.cpp)

register_test(test-adjacency-array-graph FILES test-adjacency-array-graph.cpp)

register_test(test-load-from-tree-sequence FILES test-load-from-tree-sequence.cpp LIBRARIES tskit)

register_test(test-xxhash FILES test-xxhash.cpp)

register_test(test-genomic-sequence-storage FILES test-genomic-sequence-storage.cpp LIBRARIES tskit)

register_test(
    test-allele-frequency-spectrum
    FILES
    test-allele-frequency-spectrum.cpp
    tskit-testlib/testlib.cpp
    LIBRARIES
    tskit
    fmt::fmt
)

register_test(test-allele-frequencies FILES test-allele-frequencies.cpp tskit-testlib/testlib.cpp LIBRARIES tskit)

register_test(test-compressed-forest FILES test-compressed-forest.cpp tskit-testlib/testlib.cpp LIBRARIES tskit)

register_test(
    test-tskit-site-to-tree-mapper FILES test-tskit-site-to-tree-mapper.cpp tskit-testlib/testlib.cpp LIBRARIES tskit
)

register_test(test-compressed-forest-serialization FILES test-compressed-forest-serialization.cpp LIBRARIES tskit)

register_test(test-diversity FILES test-diversity.cpp tskit-testlib/testlib.cpp LIBRARIES tskit)

register_test(test-segregating-sites FILES test-segregating-sites.cpp tskit-testlib/testlib.cpp LIBRARIES tskit)

register_test(test-tajimas-d FILES test-tajimas-d.cpp LIBRARIES tskit)

register_test(test-num-samples-below FILES test-num-samples-below.cpp)
