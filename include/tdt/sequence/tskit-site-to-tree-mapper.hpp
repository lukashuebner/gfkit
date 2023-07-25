#pragma once

#include <map>

#include <kassert/kassert.hpp>

#include "tdt/assertion_levels.hpp"
#include "tdt/graph/common.hpp"
#include "tdt/sequence/sequence.hpp"
#include "tdt/tskit.hpp"

// As the mutations are sorted by site, we can use a mapper wich does not rely on a search tree but only
// checks of the next site passes the next breakpoint in the sorted list of breakpoints.
class TSKitSiteToTreeMapper {
public:
    struct Breakpoint {
        Breakpoint(SiteId _site_id, TreeId _tree_id) noexcept : site_id(_site_id), tree_id(_tree_id) {}

        SiteId site_id;
        TreeId tree_id;
    };

    TSKitSiteToTreeMapper(TSKitTreeSequence const& tree_sequence) {
        tsk_id_t const num_sites       = asserting_cast<tsk_id_t>(tree_sequence.num_sites());
        size_t const   num_breakpoints = tree_sequence.breakpoints().size();

        _breakpoints.reserve(num_breakpoints + 1);
        SiteId site = 0;
        TreeId tree = 0;
        _breakpoints.emplace_back(site, tree);

        for (auto breakpoint: tree_sequence.breakpoints()) {
            while (site < num_sites && tree_sequence.position_of(site) < breakpoint) {
                ++site;
            }
            _breakpoints.emplace_back(site, tree);
            ++tree;
        }

        _current_breakpoint = _breakpoints.begin();
        _next_breakpoint    = _current_breakpoint + 1;
    }

    TreeId tree_id(tsk_id_t site_id) const {
        KASSERT(
            _current_breakpoint->site_id <= site_id,
            "Are the site ids not requested in ascending order?",
            tdt::assert::light
        );

        while (_next_breakpoint != _breakpoints.end() && _next_breakpoint->site_id <= site_id) {
            ++_next_breakpoint;
            ++_current_breakpoint;
        }

        return _current_breakpoint->tree_id;
    }

    TreeId operator()(tsk_id_t position) const {
        return tree_id(position);
    }

private:
    std::vector<Breakpoint>                        _breakpoints;
    mutable decltype(_breakpoints)::const_iterator _current_breakpoint;
    mutable decltype(_breakpoints)::const_iterator _next_breakpoint;
};
