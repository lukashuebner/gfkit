/*
 * MIT License
 *
 * Copyright (c) 2019-2022 Tskit Developers
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

#include "testlib.hpp"

#include <cstring>

#include <kassert/kassert.hpp>

#include "tdt/checking_casts.hpp"

char* _tmp_file_name;
FILE* _devnull;

/* Simple single tree example. */
char const* single_tree_ex_nodes = /*          6          */
    "1  0   -1   -1\n"             /*         / \         */
    "1  0   -1   -1\n"             /*        /   \        */
    "1  0   -1   -1\n"             /*       /     \       */
    "1  0   -1   -1\n"             /*      /       5      */
    "0  1   -1   -1\n"             /*     4       / \     */
    "0  2   -1   -1\n"             /*    / \     /   \    */
    "0  3   -1   -1\n";            /*   0   1   2     3   */
char const* single_tree_ex_edges = "0  1   4   0,1\n"
                                   "0  1   5   2,3\n"
                                   "0  1   6   4,5\n";
char const* single_tree_ex_sites = "0.125  0\n"
                                   "0.25   0\n"
                                   "0.5    0\n";
/* site, node, derived_state, [parent, time] */
char const* single_tree_ex_mutations = "0    2     1   -1\n"
                                       "1    4     1   -1\n"
                                       "1    0     0   1\n"  /* Back mutation over 0 */
                                       "2    0     1   -1\n" /* recurrent mutations over samples */
                                       "2    1     1   -1\n"
                                       "2    2     1   -1\n"
                                       "2    3     1   -1\n";

/*** Example from the PLOS paper ***/
/*
0.25┊     8   ┊         ┊         ┊
    ┊   ┏━┻━┓ ┊         ┊         ┊
0.20┊   ┃   ┃ ┊         ┊   7     ┊
    ┊   ┃   ┃ ┊         ┊ ┏━┻━┓   ┊
0.17┊   6   ┃ ┊   6     ┊ ┃   ┃   ┊
    ┊ ┏━┻┓  ┃ ┊ ┏━┻━┓   ┊ ┃   ┃   ┊
0.09┊ ┃  5  ┃ ┊ ┃   5   ┊ ┃   5   ┊
    ┊ ┃ ┏┻┓ ┃ ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊
0.07┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊
    ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊
0.00┊ 0 1 3 2 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊
  0.00      2.00      7.00      10.00
*/
char const* paper_ex_nodes = "1  0       -1   0\n"
                             "1  0       -1   0\n"
                             "1  0       -1   1\n"
                             "1  0       -1   1\n"
                             "0  0.071   -1   -1\n"
                             "0  0.090   -1   -1\n"
                             "0  0.170   -1   -1\n"
                             "0  0.202   -1   -1\n"
                             "0  0.253   -1   -1\n";
char const* paper_ex_edges = "2 10 4 2\n"
                             "2 10 4 3\n"
                             "0 10 5 1\n"
                             "0 2  5 3\n"
                             "2 10 5 4\n"
                             "0 7  6 0,5\n"
                             "7 10 7 0,5\n"
                             "0 2  8 2,6\n";
/* We make one mutation for each tree */
char const* paper_ex_sites     = "1      0\n"
                                 "4.5    0\n"
                                 "8.5    0\n";
char const* paper_ex_mutations = "0      2   1\n"
                                 "1      0   1\n"
                                 "2      5   1\n";
/* Two (diploid) individuals */
char const* paper_ex_individuals = "0      0.2,1.5    -1,-1\n"
                                   "0      0.0,0.0    -1,-1\n";

// Same as above but with back and recurring mutations
/*
0.25┊     8   ┊         ┊         ┊
    ┊   ┏━┻━┓ ┊         ┊         ┊
0.20┊   ┃   ┃ ┊         ┊   7     ┊
    ┊   ┃   ┃ ┊         ┊ ┏━┻━┓   ┊
0.17┊   6   ┃ ┊   6     ┊ ┃   ┃   ┊
    ┊ ┏━┻┓  ┃ ┊ ┏━┻━┓   ┊ ┃   ┃   ┊
0.09┊ ┃  5  ┃ ┊ ┃   5   ┊ ┃   5   ┊
    ┊ ┃ ┏┻┓ ┃ ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊
0.07┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊
    ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊
0.00┊ 0 1 3 2 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊
  0.00      2.00      7.00      10.00
*/
char const* multi_tree_back_recurrent_nodes = "1  0       -1   0\n"
                                              "1  0       -1   0\n"
                                              "1  0       -1   1\n"
                                              "1  0       -1   1\n"
                                              "0  0.071   -1   -1\n"
                                              "0  0.090   -1   -1\n"
                                              "0  0.170   -1   -1\n"
                                              "0  0.202   -1   -1\n"
                                              "0  0.253   -1   -1\n";

char const* multi_tree_back_recurrent_edges = "2 10 4 2\n"
                                              "2 10 4 3\n"
                                              "0 10 5 1\n"
                                              "0 2  5 3\n"
                                              "2 10 5 4\n"
                                              "0 7  6 0,5\n"
                                              "7 10 7 0,5\n"
                                              "0 2  8 2,6\n";

char const* multi_tree_back_recurrent_sites = "1      0\n"
                                              "4.5    0\n"
                                              "8.5    0\n";

/* site, node, derived_state, [parent, time] */
char const* multi_tree_back_recurrent_mutations = "0    6   1  -1\n"  // to derived state
                                                  "0    5   0   0\n"  // back to ancestral state
                                                  "0    3   1   1\n"  // and once more to the derived state
                                                  "1    5   1  -1\n"  // to derived state
                                                  "1    4   0   3\n"  // back to ancestral state
                                                  "2    4   1  -1\n"; // to derived state

/* Two (diploid) individuals */
char const* multi_tree_back_recurrent_individuals = "0      0.2,1.5    -1,-1\n"
                                               "0      0.0,0.0    -1,-1\n";


// Same as above but with back and recurring mutations and multiple derived states
/*
0.25┊     8   ┊         ┊         ┊
    ┊   ┏━┻━┓ ┊         ┊         ┊
0.20┊   ┃   ┃ ┊         ┊   7     ┊
    ┊   ┃   ┃ ┊         ┊ ┏━┻━┓   ┊
0.17┊   6   ┃ ┊   6     ┊ ┃   ┃   ┊
    ┊ ┏━┻┓  ┃ ┊ ┏━┻━┓   ┊ ┃   ┃   ┊
0.09┊ ┃  5  ┃ ┊ ┃   5   ┊ ┃   5   ┊
    ┊ ┃ ┏┻┓ ┃ ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊
0.07┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊
    ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊
0.00┊ 0 1 3 2 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊
  0.00      2.00      7.00      10.00
*/
char const* multi_derived_states_nodes = "1  0       -1   0\n"
                                         "1  0       -1   0\n"
                                         "1  0       -1   1\n"
                                         "1  0       -1   1\n"
                                         "0  0.071   -1   -1\n"
                                         "0  0.090   -1   -1\n"
                                         "0  0.170   -1   -1\n"
                                         "0  0.202   -1   -1\n"
                                         "0  0.253   -1   -1\n";

char const* multi_derived_states_edges = "2  10  4  2\n"
                                         "2  10  4  3\n"
                                         "0  10  5  1\n"
                                         "0  2   5  3\n"
                                         "2  10  5  4\n"
                                         "0  7   6  0,5\n"
                                         "7  10  7  0,5\n"
                                         "0  2   8  2,6\n";

char const* multi_derived_states_sites = "1      0\n"
                                         "4.5    0\n"
                                         "8.5    0\n";

/* site, node, derived_state, [parent, time] */
char const* multi_derived_states_mutations = "0    6   1  -1\n"  // to derived state
                                             "0    5   0   0\n"  // back to ancestral state
                                             "0    3   2   1\n"  // and once more to a derived state
                                             "1    5   3  -1\n"  // to derived state
                                             "1    4   0   3\n"  // back to ancestral state
                                             "2    4   1  -1\n"; // to derived state
// - site 0 - ancestral: 0
// node:  0  1  2  3
// state: 1  0  0  2
// - site 1 - ancestral: 0
// node:  0  1  2  3
// state: 0  3  0  0
// - site 2 - ancestral: 0
// node:  0  1  2  3
// state: 0  0  1  1

/* Two (diploid) individuals */
char const* multi_derived_states_individuals = "0      0.2,1.5    -1,-1\n"
                                               "0      0.0,0.0    -1,-1\n";

// Simple example with multiple derived states
/* Simple single tree example. */
char const* single_tree_multi_derived_states_nodes = 
                                   /*          6          */
    "1  0   -1   -1\n"             /*         / \         */
    "1  0   -1   -1\n"             /*        /   \        */
    "1  0   -1   -1\n"             /*       /     \       */
    "1  0   -1   -1\n"             /*      /       5      */
    "0  1   -1   -1\n"             /*     4       / \     */
    "0  2   -1   -1\n"             /*    / \     /   \    */
    "0  3   -1   -1\n";            /*   0   1   2     3   */
    
char const* single_tree_multi_derived_states_edges =
    "0  1   4   0,1\n"
    "0  1   5   2,3\n"
    "0  1   6   4,5\n";

char const* single_tree_multi_derived_states_sites = "0.5    0\n";

/* site, node, derived_state, [parent, time] */
char const* single_tree_multi_derived_states_mutations =
    "0    1     1   -1\n"  // To derived state 1
    "0    5     2   -1\n"; // To derived state 2
// node 0: ancestral state 0
// node 1: derived state 1
// node 2: derived state 2
// node 3: derived state 2

/* Simple utilities to parse text so we can write declarative
 * tests. This is not intended as a robust general input mechanism.
 */

void parse_nodes(char const* text, tsk_node_table_t* node_table) {
    tsk_id_t         ret_id;
    size_t           c, k;
    constexpr size_t MAX_LINE = 1024;
    char             line[MAX_LINE];
    char const*      whitespace = " \t";
    char*            p;
    double           time;
    int              population, individual;
    tsk_flags_t      flags;
    char const*      name;

    c = 0;
    while (text[c] != '\0') {
        /* Fill in the line */
        k = 0;
        while (text[c] != '\n' && text[c] != '\0') {
            KASSERT(k < MAX_LINE - 1);
            line[k] = text[c];
            c++;
            k++;
        }
        if (text[c] == '\n') {
            c++;
        }
        line[k] = '\0';
        p       = strtok(line, whitespace);
        KASSERT(p != nullptr);
        flags = asserting_cast<tsk_flags_t>(atoi(p));
        p     = strtok(NULL, whitespace);
        KASSERT(p != nullptr);
        time = atof(p);
        p    = strtok(NULL, whitespace);
        KASSERT(p != nullptr);
        population = atoi(p);
        p          = strtok(NULL, whitespace);
        if (p == nullptr) {
            individual = -1;
        } else {
            individual = atoi(p);
            p          = strtok(NULL, whitespace);
        }
        if (p == nullptr) {
            name = "";
        } else {
            name = p;
        }
        ret_id = tsk_node_table_add_row(node_table, flags, time, population, individual, name, strlen(name));
        KASSERT(ret_id >= 0);
    }
}

void parse_edges(char const* text, tsk_edge_table_t* edge_table) {
    tsk_id_t         ret_id;
    size_t           c, k;
    constexpr size_t MAX_LINE = 1024;
    char             line[MAX_LINE], sub_line[MAX_LINE];
    char const*      whitespace = " \t";
    char *           p, *q;
    double           left, right;
    tsk_id_t         parent, child;
    uint32_t         num_children;

    c = 0;
    while (text[c] != '\0') {
        /* Fill in the line */
        k = 0;
        while (text[c] != '\n' && text[c] != '\0') {
            KASSERT(k < MAX_LINE - 1);
            line[k] = text[c];
            c++;
            k++;
        }
        if (text[c] == '\n') {
            c++;
        }
        line[k] = '\0';
        p       = strtok(line, whitespace);
        KASSERT(p != nullptr);
        left = atof(p);
        p    = strtok(NULL, whitespace);
        KASSERT(p != nullptr);
        right = atof(p);
        p     = strtok(NULL, whitespace);
        KASSERT(p != nullptr);
        parent       = atoi(p);
        num_children = 0;
        p            = strtok(NULL, whitespace);
        KASSERT(p != nullptr);

        num_children = 1;
        q            = p;
        while (*q != '\0') {
            if (*q == ',') {
                num_children++;
            }
            q++;
        }
        KASSERT(num_children >= 1u);
        strncpy(sub_line, p, MAX_LINE);
        q = strtok(sub_line, ",");
        for (k = 0; k < num_children; k++) {
            KASSERT(q != nullptr);
            child  = atoi(q);
            ret_id = tsk_edge_table_add_row(edge_table, left, right, parent, child, NULL, 0);
            KASSERT(ret_id >= 0);
            q = strtok(NULL, ",");
        }
        KASSERT(q == nullptr);
    }
}

void parse_migrations(char const* text, tsk_migration_table_t* migration_table) {
    tsk_id_t         ret_id;
    size_t           c, k;
    constexpr size_t MAX_LINE = 1024;
    char             line[MAX_LINE];
    char const*      whitespace = " \t";
    char*            p;
    double           left, right, time;
    int              node, source, dest;
    char const*      metadata;

    c = 0;
    while (text[c] != '\0') {
        /* Fill in the line */
        k = 0;
        while (text[c] != '\n' && text[c] != '\0') {
            KASSERT(k < MAX_LINE - 1);
            line[k] = text[c];
            c++;
            k++;
        }
        if (text[c] == '\n') {
            c++;
        }
        line[k] = '\0';
        p       = strtok(line, whitespace);
        KASSERT(p != nullptr);
        left = atof(p);
        p    = strtok(NULL, whitespace);
        KASSERT(p != nullptr);
        right = atof(p);
        p     = strtok(NULL, whitespace);
        KASSERT(p != nullptr);
        node = atoi(p);
        p    = strtok(NULL, whitespace);
        KASSERT(p != nullptr);
        source = atoi(p);
        p      = strtok(NULL, whitespace);
        KASSERT(p != nullptr);
        dest = atoi(p);
        p    = strtok(NULL, whitespace);
        KASSERT(p != nullptr);
        time = atof(p);
        p    = strtok(NULL, whitespace);
        if (p == nullptr) {
            metadata = "";
        } else {
            metadata = p;
        }
        ret_id = tsk_migration_table_add_row(
            migration_table,
            left,
            right,
            node,
            source,
            dest,
            time,
            metadata,
            strlen(metadata)
        );
        KASSERT(ret_id >= 0);
    }
}

void parse_sites(char const* text, tsk_site_table_t* site_table) {
    tsk_id_t         ret_id;
    size_t           c, k;
    constexpr size_t MAX_LINE = 1024;
    char             line[MAX_LINE];
    double           position;
    char             ancestral_state[MAX_LINE];
    char const*      whitespace = " \t";
    char*            p;

    c = 0;
    while (text[c] != '\0') {
        /* Fill in the line */
        k = 0;
        while (text[c] != '\n' && text[c] != '\0') {
            KASSERT(k < MAX_LINE - 1);
            line[k] = text[c];
            c++;
            k++;
        }
        if (text[c] == '\n') {
            c++;
        }
        line[k] = '\0';
        p       = strtok(line, whitespace);
        KASSERT(p != nullptr);
        position = atof(p);
        p        = strtok(NULL, whitespace);
        KASSERT(p != nullptr);
        strncpy(ancestral_state, p, MAX_LINE);
        ret_id = tsk_site_table_add_row(site_table, position, ancestral_state, strlen(ancestral_state), NULL, 0);
        KASSERT(ret_id >= 0);
    }
}

void parse_mutations(char const* text, tsk_mutation_table_t* mutation_table) {
    tsk_id_t         ret_id;
    size_t           c, k;
    constexpr size_t MAX_LINE = 1024;
    char             line[MAX_LINE];
    char const*      whitespace = " \t";
    char*            p;
    tsk_id_t         node, site, parent;
    double           time;
    char             derived_state[MAX_LINE];

    /* site, node, derived_state, [parent, time] */
    c = 0;
    while (text[c] != '\0') {
        /* Fill in the line */
        k = 0;
        while (text[c] != '\n' && text[c] != '\0') {
            KASSERT(k < MAX_LINE - 1);
            line[k] = text[c];
            c++;
            k++;
        }
        if (text[c] == '\n') {
            c++;
        }
        line[k] = '\0';
        p       = strtok(line, whitespace);
        site    = atoi(p);
        KASSERT(p != nullptr);
        p = strtok(NULL, whitespace);
        KASSERT(p != nullptr);
        node = atoi(p);
        p    = strtok(NULL, whitespace);
        KASSERT(p != nullptr);
        strncpy(derived_state, p, MAX_LINE);
        parent = TSK_NULL;
        p      = strtok(NULL, whitespace);
        if (p != nullptr) {
            parent = atoi(p);
        }
        time = TSK_UNKNOWN_TIME;
        p    = strtok(NULL, whitespace);
        if (p != nullptr) {
            time = atof(p);
        }
        ret_id = tsk_mutation_table_add_row(
            mutation_table,
            site,
            node,
            parent,
            time,
            derived_state,
            strlen(derived_state),
            NULL,
            0
        );
        KASSERT(ret_id >= 0);
    }
}

void parse_individuals(char const* text, tsk_individual_table_t* individual_table) {
    tsk_id_t         ret_id;
    size_t           c, k;
    constexpr size_t MAX_LINE = 1024;
    char             line[MAX_LINE];
    char             sub_line[MAX_LINE];
    char const*      whitespace = " \t";
    char *           p, *q;
    char *           p_cont, *q_cont; // re-entrant pointers for strtok_r
    double           location[MAX_LINE];
    int              location_len;
    tsk_id_t         parents[MAX_LINE];
    int              parents_len;
    tsk_flags_t      flags;
    char const*      name;

    c = 0;
    while (text[c] != '\0') {
        /* Fill in the line */
        k = 0;
        while (text[c] != '\n' && text[c] != '\0') {
            KASSERT(k < MAX_LINE - 1);
            line[k] = text[c];
            c++;
            k++;
        }
        if (text[c] == '\n') {
            c++;
        }
        line[k] = '\0';
        p       = strtok_r(line, whitespace, &p_cont);
        KASSERT(p != nullptr);
        flags = asserting_cast<tsk_flags_t>(atoi(p));

        p = strtok_r(NULL, whitespace, &p_cont);
        KASSERT(p != nullptr);
        // the locations are comma-separated
        location_len = 1;
        q            = p;
        while (*q != '\0') {
            if (*q == ',') {
                location_len++;
            }
            q++;
        }
        KASSERT(location_len >= 1);
        strncpy(sub_line, p, MAX_LINE);
        q = strtok_r(sub_line, ",", &q_cont);
        for (k = 0; k < asserting_cast<tsk_size_t>(location_len); k++) {
            KASSERT(q != nullptr);
            location[k] = atof(q);
            q           = strtok_r(NULL, ",", &q_cont);
        }
        KASSERT(q == nullptr);

        /* parents and name are optional */
        p           = strtok_r(NULL, whitespace, &p_cont);
        parents_len = 0;
        name        = "";
        if (p != nullptr) {
            // the parents are comma-separated
            parents_len = 1;
            q           = p;
            while (*q != '\0') {
                if (*q == ',') {
                    parents_len++;
                }
                q++;
            }
            KASSERT(parents_len >= 1);
            strncpy(sub_line, p, MAX_LINE);
            q = strtok_r(sub_line, ",", &q_cont);
            for (k = 0; k < asserting_cast<tsk_size_t>(parents_len); k++) {
                KASSERT(q != nullptr);
                parents[k] = atoi(q);
                q          = strtok_r(NULL, ",", &q_cont);
            }
            KASSERT(q == nullptr);
            p = strtok_r(NULL, whitespace, &p_cont);
            if (p != nullptr) {
                name = p;
            }
        }
        ret_id = tsk_individual_table_add_row(
            individual_table,
            flags,
            location,
            asserting_cast<tsk_size_t>(location_len),
            parents,
            asserting_cast<tsk_size_t>(parents_len),
            name,
            strlen(name)
        );
        KASSERT(ret_id >= 0);
    }
}

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
) {
    int                    ret;
    tsk_id_t               ret_id;
    tsk_table_collection_t tables;
    tsk_id_t               max_population_id;
    tsk_size_t             j;

    KASSERT(ts != nullptr);
    KASSERT(nodes != nullptr);
    KASSERT(edges != nullptr);
    /* Not supporting provenance here for now */
    KASSERT(provenance == nullptr);

    ret = tsk_table_collection_init(&tables, tc_options);
    KASSERT(ret == 0);
    tables.sequence_length = sequence_length;
    parse_nodes(nodes, &tables.nodes);
    parse_edges(edges, &tables.edges);
    if (sites != nullptr) {
        parse_sites(sites, &tables.sites);
    }
    if (mutations != nullptr) {
        parse_mutations(mutations, &tables.mutations);
    }
    if (individuals != nullptr) {
        parse_individuals(individuals, &tables.individuals);
    }
    if (migrations != nullptr) {
        parse_migrations(migrations, &tables.migrations);
    }
    /* We need to add in populations if they are referenced */
    max_population_id = -1;
    for (j = 0; j < tables.nodes.num_rows; j++) {
        max_population_id = TSK_MAX(max_population_id, tables.nodes.population[j]);
    }
    for (j = 0; j < tables.migrations.num_rows; j++) {
        max_population_id = TSK_MAX(max_population_id, tables.migrations.source[j]);
        max_population_id = TSK_MAX(max_population_id, tables.migrations.dest[j]);
    }
    if (max_population_id >= 0) {
        for (j = 0; j <= (tsk_size_t)max_population_id; j++) {
            ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
            KASSERT(asserting_cast<decltype(j)>(ret_id) == j);
        }
    }

    ret = tsk_table_collection_sort(&tables, 0, 0);
    KASSERT(ret == 0);
    ret = tsk_table_collection_build_index(&tables, 0);
    KASSERT(ret == 0);
    if (mutations != NULL) {
        ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
        KASSERT(ret == 0);
    }
    ret = tsk_treeseq_init(ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    /* tsk_treeseq_print_state(ts, stdout); */
    /* printf("ret = %s\n", tsk_strerror(ret)); */
    KASSERT(ret == 0);
    tsk_table_collection_free(&tables);
}
