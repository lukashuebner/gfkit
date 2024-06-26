#pragma once

#include <kassert/kassert.hpp>
#include <tskit/core.h>

#include "sfkit/assertion_levels.hpp"
#include "sfkit/graph/TsToSfNodeMapper.hpp"
#include "sfkit/sequence/GenomicSequence.hpp"
#include "sfkit/tskit/tskit.hpp"

namespace sfkit::sequence {
// We interleaf building the genomic sequence storage with building the forest to save memory.
// To build the genomic sequence storage, we need the mapping of ts nodes to sf subtree ids for all trees.
// If we were to build the forest first, we would need to store this mapping for all trees, which takes up hundreds
// of gigabytes of memory. We therefore iterate over the trees only once, first extending the forest, then
// processing the mutations for the current tree and only then advancing to the next tree in the tree sequence.
class GenomicSequenceFactory {
public:
    // Rename storage to store?
    GenomicSequenceFactory(tskit::TSKitTreeSequence const& tree_sequence)
        : _sequence(tree_sequence.num_sites(), tree_sequence.num_mutations()),
          _site2tree(tree_sequence),
          _mutation_it(tree_sequence.mutations().begin()),
          _mutations_end(tree_sequence.mutations().end()) {
        _set_ancestral_states(tree_sequence);
    }

    // Call this for all trees in order, make sure the the mutations are sorted by site.
    // return true if done; else returns false
    template <typename TsToSfNodeMapper>
    bool process_mutations(TreeId tree_id, TsToSfNodeMapper const& ts_to_sf_node) {
        // KASSERT(
        //     ts_node2sf_subtree.num_requested() == ts_node2sf_subtree.size(),
        //     "The number of requested and mapped ts node IDs differ.",
        //     sfkit::assert::light
        // );
        while (_mutation_it != _mutations_end) {
            // TODO Use less bits for site and tree ids
            tsk_id_t const site_id       = _mutation_it->site;
            TreeId const   sites_tree_id = _site2tree(site_id);
            KASSERT(
                sites_tree_id >= tree_id,
                "We seemed to have missed processing a mutation. Are the mutations sorted by tree id?",
                sfkit::assert::light
            );
            if (sites_tree_id > tree_id) {
                return false; // Continue later, first the Forest compressor needs to extend the DAG and compute the new
                              // mapping of ts nodes to sf subtree ids.
            }

            NodeId const sf_node_id = ts_to_sf_node(asserting_cast<size_t>(_mutation_it->node));
            KASSERT(_mutation_it->node != TSK_NULL, "Mutation node is null", sfkit::assert::light);
            KASSERT(_mutation_it->derived_state_length == 1u, "Derived state length is not 1", sfkit::assert::light);

            AllelicState const derived_state = *_mutation_it->derived_state;
            tsk_id_t const     mutation_id   = asserting_cast<tsk_id_t>(_sequence.num_mutations());
            KASSERT(
                mutation_id == _mutation_it->id,
                "Mutation ID is not equal to the index in the mutations vector",
                sfkit::assert::light
            );

            tsk_id_t const     parent_mutation_id = _mutation_it->parent;
            AllelicState const ancestral_state =
                parent_mutation_id == TSK_NULL
                    ? _sequence.ancestral_state(asserting_cast<SiteId>(site_id))
                    : _sequence.mutation_by_id(asserting_cast<MutationId>(parent_mutation_id)).allelic_state();

            _sequence.emplace_back(site_id, tree_id, sf_node_id, derived_state, ancestral_state);
            ++_mutation_it;
        }
        return true;
    }

    GenomicSequence&& move_storage() {
        KASSERT(
            _finalized,
            "Storage has not been finalized yet (we need to build the mutation indiced).",
            sfkit::assert::light
        );
        KASSERT(!_moved, "Storage has already been moved", sfkit::assert::light);
        _moved = true;
        return std::move(_sequence);
    }

    void finalize() {
        KASSERT(!_finalized, "Storage has already been finalized", sfkit::assert::light);
        _finalized = true;
        _sequence.build_mutation_indices();
    }

private:
    GenomicSequence                  _sequence;
    tskit::TskMutationView           _tsk_mutations;
    TSKitSiteToTreeMapper            _site2tree;
    tskit::TskMutationView::iterator _mutation_it;
    tskit::TskMutationView::iterator _mutations_end;
    bool                             _finalized = false;
    bool                             _moved     = false;

    void _set_ancestral_states(tskit::TSKitTreeSequence const& tree_sequence) {
        // Store ancestral states
        for (auto&& site: tree_sequence.sites()) {
            KASSERT(site.ancestral_state_length == 1u, "Ancestral state length is not 1", sfkit::assert::light);
            _sequence.emplace_back(*site.ancestral_state);
        }
        KASSERT(
            _sequence.num_sites() == tree_sequence.num_sites(),
            "Number of sites reported by num_sites() and in the sites() iterator does not match",
            sfkit::assert::light
        );
    }
};
} // namespace sfkit::sequence
