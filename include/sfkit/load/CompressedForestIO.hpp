#include <fstream>

#include <fmt/format.h>

#if defined(__GNUC__) && !defined(__clang__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wnoexcept"
#endif
#include <cereal/archives/binary.hpp>
#if defined(__GNUC__) && !defined(__clang__)
    #pragma GCC diagnostic pop
#endif

#include "sfkit/SuccinctForest.hpp"
#include "sfkit/graph/CompressedForest.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"

class CompressedForestIO {
public:
    using Version = uint64_t;

    static void load(std::string const& filename, CompressedForest& forest, GenomicSequence& sequence) {
        std::ifstream              is(filename, std::ios::binary | std::ios::in);
        cereal::BinaryInputArchive archive(is);

        Version archive_version;
        archive(archive_version, forest, sequence);

        if (archive_version != CURRENT_VERSION) {
            throw std::runtime_error(fmt::format(
                "Archive {} has version {} but current version is {}",
                filename,
                archive_version,
                CURRENT_VERSION
            ));
        }

        KASSERT(
            forest.num_nodes_is_set(),
            "Number of nodes of the loaded forest are not computed",
            sfkit::assert::light
        );
        KASSERT(
            sequence.mutation_indices_are_built(),
            "Mutation indices of the loaded sequence are not built",
            sfkit::assert::light
        );
    }

    static void save(std::string const& filename, CompressedForest& forest, GenomicSequence& sequence) {
        std::ofstream               os(filename, std::ios::binary | std::ios::out);
        cereal::BinaryOutputArchive archive(os);

        // Build mutation indices and compute the nodes of the (edge list) DAG
        // The mutation indices are build during the serialization of the sequence
        // The number of nodes of the DAG are computed during the deserialization of the forest
        archive(CURRENT_VERSION, forest, sequence);
    }

    // static SequenceForest load(std::string const& filename, SequenceForest& sequence_forest) {
    //     SequenceForest sequence_forest();
    //     load(filename, sequence_forest._forest, sequence_forest._sequence);
    // }

    // static void save(std::string& filename, SequenceForest& sequence_forest) {
    //     save(filename, sequence_forest._forest, sequence_forest._sequence);
    // }

private:
    static constexpr Version CURRENT_VERSION = 2;
};
