#pragma once

template <typename T>
concept CompressedForestC = requires(T t) {
    { t.num_nodes() } -> std::convertible_to<NodeId>;
    { all_samples() } -> std::convertible_to<SampleSet>;

SampleId num_samples() const {
TreeId num_trees() const {
SampleId num_leaves() const {
NodeId num_unique_subtrees() const {
NodeId num_backrefs() const {
bool is_leaf(size_t const bp_idx) const {
SampleId leaf_idx_to_id(SampleId const leaf_id) const {
NodeId node_id(size_t const bp_idx) const {
[[nodiscard]] bool operator==(BPCompressedForest const& other) {

private:
// TODO Which of these vectors are actually needed?
// TODO Use compressed bit-vectors
sdsl::bit_vector        _is_reference;
sdsl::rank_support_v5<> _is_reference_rank;
sdsl::bit_vector        _is_leaf;
sdsl::rank_support_v5<> _is_leaf_rank;
sdsl::bit_vector        _balanced_parenthesis;
// TODO Document why we're ranking the 0 pattern
sdsl::rank_support_v5<0>          _balanced_parenthesis_rank;
sdsl::int_vector<NodeId_bitwidth> _references;
sdsl::int_vector<NodeId_bitwidth> _leaves;
NodeId                            _num_nodes;
SampleId                          _num_leaves;
TreeId                            _num_trees;
};

} // namespace sfkit::bp

} // namespace sfkit::dag

};
