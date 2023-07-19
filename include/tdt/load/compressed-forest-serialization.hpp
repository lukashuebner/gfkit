#include "fmt/format.h"
#include <fstream>

#include <cereal/archives/binary.hpp>

#include "tdt/graph/compressed-forest.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"

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

        KASSERT(forest.nodes_are_computed(), "Nodes of the loaded forest are not computed", tdt::assert::light);
        KASSERT(
            sequence.mutation_indices_are_built(),
            "Mutation indices of the loaded sequence are not built",
            tdt::assert::light
        );
    }

    static void save(std::string const& filename, CompressedForest& forest, GenomicSequence& sequence) {
        std::ofstream               os(filename, std::ios::binary | std::ios::out);
        cereal::BinaryOutputArchive archive(os);

        // Build mutation indices and compute the nodes of the (edge list) DAG
        // The mutation indices are build during the serialization of the sequence
        // The nodes of the DAG are computed during the deserialization of the forest
        archive(CURRENT_VERSION, forest, sequence);
    }

private:
    static constexpr Version CURRENT_VERSION = 1;
};
