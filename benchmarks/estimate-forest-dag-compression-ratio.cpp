#define ENABLE_TIMERS

#include <cassert>
#include <filesystem>
#include <iostream>
#include <string>

#include <openssl/evp.h>

#include "timer.hpp"

using namespace genesis::tree;
using namespace genesis::utils;

class SubtreeId {
public:
    SubtreeId() = default;

    SubtreeId(std::string& leaf_label) {
        Sha256 digester;
        digester.update(leaf_label);
        _subtree_digest = digester.digest();
        assert(_subtree_digest.size() == 32);
    }

    SubtreeId(std::vector<SubtreeId>& children) {
        Sha256 digester;
        for (auto& child: children) {
            digester.update(child._subtree_digest);
        }
        _subtree_digest = digester.digest();
        assert(_subtree_digest.size() == 32);
    }

    bool operator==(SubtreeId const& other) const {
        return _subtree_digest == other._subtree_digest;
    }

    bool operator!=(SubtreeId const& other) const {
        return !(*this == other);
    }

    Sha256Digest digest() const {
        assert(_subtree_digest.size() == 32);
        return _subtree_digest;
    }

    std::string to_string() const {
        std::string digest_str;
        for (auto& byte: _subtree_digest) {
            digest_str += std::to_string(byte);
        }
        return digest_str;
    }

private:
    Sha256Digest _subtree_digest;
};


using VertexId = size_t;

class Edge {
public:
    Edge(VertexId from, VertexId to) : _from(from), _to(to) {}
    VertexId from() const {
        return _from;
    }
    VertexId to() const {
        return _to;
    }

private:
    const VertexId _from;
    const VertexId _to;
};

using EdgeList = std::vector<Edge>;

class EdgeListGraph {
public:
    EdgeListGraph() {}

    void add_edge(VertexId from, VertexId to) {
        _edges.emplace_back(Edge(from, to));
    }
    void add_root(VertexId root) {
        _roots.push_back(root);
    }
    void add_leaf(VertexId leaf) {
        _leaves.push_back(leaf);
    }
    size_t size() const {
        return _edges.size();
    }

    // template <typename Visitor> void visit(Visitor &&visitor) const {
    //   for (auto &edge : _edges) {
    //     visitor(edge.from(), edge.to());
    //   }
    // }

    using iterator       = EdgeList::iterator;
    using const_iterator = EdgeList::const_iterator;

    iterator begin() {
        return _edges.begin();
    }
    iterator end() {
        return _edges.end();
    }
    const_iterator begin() const {
        return _edges.begin();
    }
    const_iterator end() const {
        return _edges.end();
    }
    const_iterator cbegin() const {
        return _edges.cbegin();
    }
    const_iterator cend() const {
        return _edges.cend();
    }

    std::vector<NodeId> const& roots() const {
        return _roots;
    }
    std::vector<NodeId> const& leaves() const {
        return _leaves;
    }

    bool directed() const {
        return true;
    }

    // TODO Implement roots_have_in_edges()

private:
    EdgeList            _edges;
    std::vector<NodeId> _roots;
    std::vector<NodeId> _leaves;
};

class SubtreeIdNodeMapper {
public:
    SubtreeIdNodeMapper() {}

    VertexId insert(SubtreeId subtree_id) {
        assert(_subtree_to_node_map.find(subtree_id) == _subtree_to_node_map.end());
        _subtree_to_node_map[subtree_id] = _next_node_id;
        return _next_node_id++;
    }

    // TODO Implement get with optional (if not found, return nullopt)
    bool contains(SubtreeId subtree_id) const {
        return _subtree_to_node_map.find(subtree_id) != _subtree_to_node_map.end();
    }

    VertexId get(SubtreeId subtree_id) const {
        assert(_subtree_to_node_map.find(subtree_id) != _subtree_to_node_map.end());
        return _subtree_to_node_map.at(subtree_id);
    }

    VertexId operator[](SubtreeId subtree_id) const {
        return get(subtree_id);
    }

    VertexId num_nodes() const {
        return _next_node_id;
    }

private:
    std::unordered_map<SubtreeId, NodeId> _subtree_to_node_map;
    VertexId                              _next_node_id = 0;
};

int main() {
    const std::string DATA_DIR = "../../../data/tgp_chr08/testset";

    std::vector<Tree>      trees;
    std::vector<SubtreeId> subtrees;
    for (auto const& entry: std::filesystem::directory_iterator(DATA_DIR)) {
        if (entry.path().extension() == ".newick") {
            std::cerr << "Reading " << entry.path().filename() << std::endl;
            trees.push_back(CommonTreeNewickReader().read(from_file(entry.path())));
            auto& tree = trees.back();

            std::vector<bool> node_visited(tree.node_count(), false);
            std::for_each(node_visited.begin(), node_visited.end(), [](bool visited) { assert(!visited); });

#ifndef NDEBUG
            [[maybe_unused]] size_t reachable_from_root = 0;
            for (auto it: postorder(tree)) {
                reachable_from_root++;
                node_visited[it.node().index()] = true;
            }
            assert(reachable_from_root == tree.node_count());
            std::for_each(node_visited.begin(), node_visited.end(), [](bool visited) { assert(visited); });

            size_t root_index = tree.root_node().index();
            for (auto it: postorder(tree)) {
                TreeLink* link = &it.node().primary_link();
                while (link->node().index() != root_index) {
                    link = &link->node().primary_link();
                }
            }
#endif
        }
    }

    std::unordered_set<SubtreeId> subtree_ids;
    uint64_t                      total_subtrees = 0;

    EdgeListGraph       postorder_edges;
    SubtreeIdNodeMapper subtree_to_node_map;
    size_t              tree_no = 0;

    for (auto& tree: trees) {
        std::cerr << "Processing tree " << tree_no++ << std::endl;
        std::cout << "Read: Tree has " << tree.node_count() << " nodes" << std::endl;
        std::vector<SubtreeId> this_trees_subtree_ids(tree.node_count());
        std::for_each(this_trees_subtree_ids.begin(), this_trees_subtree_ids.end(), [](SubtreeId& subtree_id) {
            assert(subtree_id == SubtreeId());
        });
        std::vector<bool> node_visited(tree.node_count(), false);
        std::for_each(node_visited.begin(), node_visited.end(), [](bool visited) { assert(!visited); });

        for (auto it: postorder(tree)) {
            size_t node_index = it.node().index();
            assert(!node_visited[node_index]);
            node_visited[node_index] = true;

            total_subtrees++;
            auto* link = &(it.node().primary_link());
            // Assert that the node on the edge towards the root has not been visited
            // yet
            assert(is_root(it.node().primary_link()) || !node_visited[it.node().primary_link().outer().node().index()]);

            if (&(link->next()) == link) { // Node is leaf
                // Compute the subtree id of this leaf node by hashing its name
                SubtreeId subtree_id(link->node().data<CommonNodeData>().name);

                // Store the subtree id for this node
                assert(link->node().index() == node_index);
                assert(this_trees_subtree_ids.size() > node_index);
                this_trees_subtree_ids[node_index] = subtree_id;

                // Add the subtree id to the set of all subtree ids
                subtree_ids.insert(subtree_id);

                // Map the subtree id to a node id to be used in the DAG
                if (!subtree_to_node_map.contains(subtree_id)) {
                    NodeId node_id = subtree_to_node_map.insert(subtree_id);
                    postorder_edges.add_leaf(node_id);
                    assert(link->outer().outer().node().index() == node_index);
                    assert(is_root(link->outer()) || &link->outer() != &link->outer().node().primary_link());
                }
            } else { // Node is inner node
                // Compute the subtree id of this inner node by hashing the subtree ids
                // of it's children
                std::vector<SubtreeId> children;
                std::vector<NodeId>    children_node_ids;
                auto*                  link_to_root = link;
                auto*                  child        = &(link->next());
                if (is_root(*link)) {
                    child = link_to_root;
                }
                do {
                    auto child_index = child->outer().node().index();
                    assert(node_visited[child_index]);
                    auto child_subtree_id = this_trees_subtree_ids[child_index];
                    children.emplace_back(child_subtree_id);
                    assert(subtree_to_node_map.contains(child_subtree_id));
                    children_node_ids.emplace_back(subtree_to_node_map[child_subtree_id]);
                    child = &(child->next());
                } while (link_to_root != child);
                SubtreeId subtree_id(children);

                // Add the subtree id to the set of all subtree ids
                subtree_ids.insert(subtree_id);

                // Map the subtree id to a node id to be used in the DAG
                assert(link->node().index() == node_index);
                assert(this_trees_subtree_ids.size() > node_index);
                this_trees_subtree_ids[node_index] = subtree_id;

                // Add this node to the DAG if not already present. As the DAG is stored
                // as a list of edges, we need to add an edge from this node to each of
                // it's children.
                if (!subtree_to_node_map.contains(subtree_id)) {
                    NodeId node_id = subtree_to_node_map.insert(subtree_id);
                    for (auto child_node_id: children_node_ids) {
                        postorder_edges.add_edge(node_id, child_node_id);
                    }
                }
            }
        }

        // Add this tree's root to the list of roots.
        auto root_node_id = subtree_to_node_map[this_trees_subtree_ids[tree.root_node().index()]];
        postorder_edges.add_root(root_node_id);
    }
    uint64_t unique_subtrees = subtree_ids.size();
    assert(unique_subtrees == subtree_to_node_map.num_nodes());
    std::cerr << "Unique subtrees (== nodes in DAG): " << unique_subtrees << std::endl;
    std::cerr << "Total number of subtrees: " << total_subtrees << std::endl;
    assert(postorder_edges.roots().size() == trees.size());
    std::cerr << "Number of trees: " << trees.size() << std::endl;
    std::cerr << "Number of edges: " << postorder_edges.size() << std::endl;

    // Compute the number of samples below each root.
    TIME_NEXT_SECTION("ComputeSubtreeSizes");
    do_not_optimize(postorder_edges);
    std::vector<size_t> subtree_sizes(subtree_to_node_map.num_nodes(), 0);
    for (size_t leaf: postorder_edges.leaves()) {
        subtree_sizes[leaf] = 1;
    }
    std::vector<bool> pointed_to(subtree_to_node_map.num_nodes(), false);
    for (auto& edge: postorder_edges) {
        assert(!pointed_to[edge.from()]);
        pointed_to[edge.to()] = true;
        subtree_sizes[edge.from()] += subtree_sizes[edge.to()];
    }
    do_not_optimize(subtree_sizes);
    TIME_STOP();
    for (auto root: postorder_edges.roots()) {
        // assert(!pointed_to[root]);
        pointed_to[root] = true;
    }
    size_t num_pointed_to = std::count(pointed_to.begin(), pointed_to.end(), true);
    std::cerr << "Number of nodes pointed to: " << num_pointed_to << std::endl;
    //   std::for_each(pointed_to.begin(), pointed_to.end(),
    //                 [](bool b) { assert(b); });
    TimerRegister::instance().print("info");

#if not(defined(NDEBUG))
    for (VertexId node_id = 0; node_id < subtree_to_node_map.num_nodes(); node_id++) {
        if (!pointed_to[node_id]) {
            std::cerr << "No one points to node " << node_id << std::endl;
        }
    }
#endif

    for (auto root: postorder_edges.roots()) {
        std::cout << "Computed: Tree has " << subtree_sizes[root] << " nodes" << std::endl;
    }

    return 0;
}
