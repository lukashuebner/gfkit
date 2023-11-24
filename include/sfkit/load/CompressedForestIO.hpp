#include <fstream>

#include <fmt/format.h>

#include "sfkit/SuccinctForest.hpp"
#include "sfkit/graph/CompressedForest.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"

namespace sfkit::io::internal {
using Version = uint64_t;
using Magic   = uint64_t;

static constexpr Version DAG_ARCHIVE_VERSION = 3;
static constexpr Magic   DAG_ARCHIVE_MAGIC   = 1307950585415129820;

static constexpr Version BP_ARCHIVE_VERSION = 1;
static constexpr Magic   BP_ARCHIVE_MAGIC   = 7612607674453629763;

class DAGCompressedForestIO {
public:
    static void load(std::string const& filename, CompressedForest& forest, GenomicSequence& sequence) {
        std::ifstream              is(filename, std::ios::binary | std::ios::in);
        cereal::BinaryInputArchive archive(is);

        Magic   magic;
        Version version;
        archive(magic, version, forest, sequence);

        if (magic != DAG_ARCHIVE_MAGIC) {
            throw std::runtime_error("Archive is invalid (mismatch of magic number)");
        }

        if (version != DAG_ARCHIVE_VERSION) {
            throw std::runtime_error(fmt::format(
                "Archive {} has version {} but current version is {}",
                filename,
                version,
                DAG_ARCHIVE_VERSION
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
        std::ofstream os(filename, std::ios::binary | std::ios::out);

        cereal::BinaryOutputArchive archive(os);
        archive(DAG_ARCHIVE_MAGIC, DAG_ARCHIVE_VERSION, forest, sequence);

        os.close();
    }
};

class BPCompressedForestIO {
public:
    static void save(std::string const& filename, BPCompressedForest& forest, GenomicSequence& sequence) {
        std::ofstream os(filename, std::ios::binary | std::ios::out);

        // TODO Document decicion not to use cereal here
        os.write(reinterpret_cast<char const*>(&BP_ARCHIVE_MAGIC), sizeof(BP_ARCHIVE_MAGIC));
        os.write(reinterpret_cast<char const*>(&BP_ARCHIVE_VERSION), sizeof(BP_ARCHIVE_VERSION));

        forest.save(os);
        sequence.save(os);
    }

    static void load(std::string const& filename, BPCompressedForest& forest, GenomicSequence& sequence) {
        std::ifstream is(filename, std::ios::binary | std::ios::in);

        Magic   magic;
        Version version;

        is.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        is.read(reinterpret_cast<char*>(&version), sizeof(version));

        if (magic != BP_ARCHIVE_MAGIC) {
            throw std::runtime_error("Archive is invalid (mismatch of magic number)");
        }

        if (version != BP_ARCHIVE_VERSION) {
            throw std::runtime_error(fmt::format(
                "Archive {} has version {} but current version is {}",
                filename,
                version,
                BP_ARCHIVE_VERSION
            ));
        }

        forest.load(is);
        sequence.load(is);

        KASSERT(
            sequence.mutation_indices_are_built(),
            "Mutation indices of the loaded sequence are not built",
            sfkit::assert::light
        );
    }
};
} // namespace sfkit::io::internal

namespace sfkit::io {
class CompressedForestIO {
public:
    static void load(std::string const& filename, BPCompressedForest& forest, GenomicSequence& sequence) {
        internal::BPCompressedForestIO::load(filename, forest, sequence);
    }

    static void save(std::string const& filename, BPCompressedForest& forest, GenomicSequence& sequence) {
        internal::BPCompressedForestIO::save(filename, forest, sequence);
    }

    static void load(std::string const& filename, CompressedForest& forest, GenomicSequence& sequence) {
        internal::DAGCompressedForestIO::load(filename, forest, sequence);
    }

    static void save(std::string const& filename, CompressedForest& forest, GenomicSequence& sequence) {
        internal::DAGCompressedForestIO::save(filename, forest, sequence);
    }
};
} // namespace sfkit::io
