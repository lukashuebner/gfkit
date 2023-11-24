/*
 * MIT License
 *
 * Copyright (c) 2019-2021 Tskit Developers
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#ifndef __TESTLIB_H__
    #define __TESTLIB_H__

    #include "tskit/core.h"
    #include "tskit/tables.h"

    #include <stdio.h>
    #include <stdlib.h>
    #include <tskit/trees.h>
    #include <unistd.h>

/* Global variables used in the test suite */

extern char* _tmp_file_name;
extern FILE* _devnull;

void tsk_treeseq_from_text(
    tsk_treeseq_t* ts,
    double         sequence_length,
    char const*    nodes,
    char const*    edges,
    char const*    migrations,
    char const*    sites,
    char const*    mutations,
    char const*    individuals,
    char const*    provenance,
    tsk_flags_t    tc_options
);

void parse_nodes(char const* text, tsk_node_table_t* node_table);
void parse_edges(char const* text, tsk_edge_table_t* edge_table);
void parse_sites(char const* text, tsk_site_table_t* site_table);
void parse_mutations(char const* text, tsk_mutation_table_t* mutation_table);
void parse_individuals(char const* text, tsk_individual_table_t* individual_table);

extern char const* single_tree_ex_nodes;
extern char const* single_tree_ex_edges;
extern char const* single_tree_ex_sites;
extern char const* single_tree_ex_mutations;

extern char const* multiple_tree_ex_nodes;
extern char const* multiple_tree_ex_edges;

extern char const* empty_ex_nodes;
extern char const* empty_ex_edges;

extern char const* paper_ex_nodes;
extern char const* paper_ex_edges;
extern char const* paper_ex_sites;
extern char const* paper_ex_mutations;
extern char const* paper_ex_individuals;

extern char const* multi_tree_back_recurrent_nodes;
extern char const* multi_tree_back_recurrent_edges;
extern char const* multi_tree_back_recurrent_sites;
extern char const* multi_tree_back_recurrent_mutations;
extern char const* multi_tree_back_recurrent_individuals;

extern char const* multi_derived_states_nodes;
extern char const* multi_derived_states_edges;
extern char const* multi_derived_states_sites;
extern char const* multi_derived_states_mutations;
extern char const* multi_derived_states_individuals;

extern char const* single_tree_multi_derived_states_nodes;
extern char const* single_tree_multi_derived_states_edges;
extern char const* single_tree_multi_derived_states_sites;
extern char const* single_tree_multi_derived_states_mutations;

#endif
