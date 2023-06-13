#include <fstream>

#include <cereal/archives/binary.hpp>

#include "tdt/graph/compressed-forest.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"

class CompressedForestIO {
public:
    static void
    load(std::string const& filename, CompressedForest& compressed_forest, GenomicSequenceStorage& sequence_store) {
        std::ifstream              is(filename, std::ios::binary | std::ios::in);
        cereal::BinaryInputArchive archive(is);

        archive(compressed_forest, sequence_store);
    }

    static void
    save(std::string const& filename, CompressedForest& compressed_forest, GenomicSequenceStorage& sequence_store) {
        std::ofstream               os(filename, std::ios::binary | std::ios::out);
        cereal::BinaryOutputArchive archive(os);

        archive(compressed_forest, sequence_store);
    }
};
