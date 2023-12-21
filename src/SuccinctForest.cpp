#include "sfkit/SuccinctForest.hpp"

#include "sfkit/bp/BPCompressedForest.hpp"
#include "sfkit/dag/DAGCompressedForest.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"

namespace sfkit {

using DAGSuccinctForest        = SuccinctForest<DAGCompressedForest, PerfectDNAHasher>;
using DAGSuccinctForestNumeric = SuccinctForest<DAGCompressedForest, PerfectNumericHasher>;
using BPSuccinctForest         = SuccinctForest<BPCompressedForest, PerfectDNAHasher>;
using BPSuccinctForestNumeric  = SuccinctForest<BPCompressedForest, PerfectNumericHasher>;

template class SuccinctForest<DAGCompressedForest, PerfectDNAHasher>;
template class SuccinctForest<DAGCompressedForest, PerfectNumericHasher>;
template class SuccinctForest<BPCompressedForest, PerfectDNAHasher>;
template class SuccinctForest<BPCompressedForest, PerfectNumericHasher>;

} // namespace sfkit
