#include <fstream>

#include <cereal/archives/binary.hpp>

#include "tdt/graph/compressed-forest.hpp"
#include "tdt/sequence/genomic-sequence-storage.hpp"

class CompressedForestIO {
public:
    static void load(std::string const& filename, CompressedForest& forest, GenomicSequence& sequence) {
        std::ifstream              is(filename, std::ios::binary | std::ios::in);
        cereal::BinaryInputArchive archive(is);

        archive(forest, sequence);
        sequence.build_mutation_indices();
    }

    static void save(std::string const& filename, CompressedForest& forest, GenomicSequence& sequence) {
        std::ofstream               os(filename, std::ios::binary | std::ios::out);
        cereal::BinaryOutputArchive archive(os);

        archive(forest, sequence);
    }
};
