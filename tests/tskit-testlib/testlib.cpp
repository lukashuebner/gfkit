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
// char const* multi_tree_back_recurrent_mutations = "0    6   1  -1\n"  // to derived state
//                                                   "0    5   0   0\n"  // back to ancestral state
//                                                   "0    3   1   1\n"  // and once more to the derived state
//                                                   "1    5   1  -1\n"  // to derived state
//                                                   "1    4   0   3\n"  // back to ancestral state
//                                                   "2    4   1  -1\n"; // to derived state

char const* multi_tree_back_recurrent_mutations = "0    6   1\n"  // to derived state
                                                  "0    5   0\n"  // back to ancestral state
                                                  "0    3   1\n"  // and once more to the derived state
                                                  "1    5   1\n"  // to derived state
                                                  "1    4   0\n"  // back to ancestral state
                                                  "2    4   1\n"; // to derived state

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

char const* multi_derived_states_edges = "2 10 4 2\n"
                                         "2 10 4 3\n"
                                         "0 10 5 1\n"
                                         "0 2  5 3\n"
                                         "2 10 5 4\n"
                                         "0 7  6 0,5\n"
                                         "7 10 7 0,5\n"
                                         "0 2  8 2,6\n";

char const* multi_derived_states_sites = "1      0\n"
                                         "4.5    0\n"
                                         "8.5    0\n";

/* site, node, derived_state, [parent, time] */
// char const* multi_tree_back_recurrent_mutations = "0    6   1  -1\n"  // to derived state
//                                                   "0    5   0   0\n"  // back to ancestral state
//                                                   "0    3   1   1\n"  // and once more to the derived state
//                                                   "1    5   1  -1\n"  // to derived state
//                                                   "1    4   0   3\n"  // back to ancestral state
//                                                   "2    4   1  -1\n"; // to derived state

char const* multi_derived_states_mutations = "0    6   1\n"  // to derived state
                                             "0    5   0\n"  // back to ancestral state
                                             "0    3   2\n"  // and once more to the derived state
                                             "1    5   3\n"  // to derived state
                                             "1    4   0\n"  // back to ancestral state
                                             "2    4   1\n"; // to derived state

/* Two (diploid) individuals */
char const* multi_derived_states_individuals = "0      0.2,1.5    -1,-1\n"
                                               "0      0.0,0.0    -1,-1\n";

/*** An example of a nonbinary tree sequence ***/
/*
0.41┊          12     ┊         12      ┊
    ┊           ┃     ┊          ┃      ┊
0.28┊           ┃     ┊          ┃      ┊
    ┊           ┃     ┊          ┃      ┊
0.13┊          10     ┊         10      ┊
    ┊         ┏━╋━┓   ┊         ┏┻┓     ┊
0.07┊         ┃ ┃ ┃   ┊         ┃ ┃     ┊
    ┊         ┃ ┃ ┃   ┊         ┃ ┃     ┊
0.01┊         ┃ ┃ ┃   ┊         ┃ ┃     ┊
    ┊         ┃ ┃ ┃   ┊         ┃ ┃     ┊
0.00┊ 0 1 2 3 4 5 7 6 ┊ 0 1 2 3 4 7 5 6 ┊
*/
char const* nonbinary_ex_nodes     = "1  0       0   -1\n"
                                     "1  0       0   -1\n"
                                     "1  0       0   -1\n"
                                     "1  0       0   -1\n"
                                     "1  0       0   -1\n"
                                     "1  0       0   -1\n"
                                     "1  0       0   -1\n"
                                     "1  0       0   -1\n"
                                     "0  0.01    0   -1\n"
                                     "0  0.068   0   -1\n"
                                     "0  0.130   0   -1\n"
                                     "0  0.279   0   -1\n"
                                     "0  0.405   0   -1\n";
char const* nonbinary_ex_edges     = "0	100	8	0,1,2,3\n"
                                     "0	100	9	6,8\n"
                                     "0  100 10  4\n"
                                     "0  17  10  5\n"
                                     "0  100 10  7\n"
                                     "17	100	11	5,9\n"
                                     "0	17	12	9\n"
                                     "0  100 12  10\n"
                                     "17	100	12	11";
char const* nonbinary_ex_sites     = "1  0\n"
                                     "18 0\n";
char const* nonbinary_ex_mutations = "0    2   1\n"
                                     "1    11  1";

/*** An example of a tree sequence with unary nodes. ***/
/*
0.25┊     8   ┊   8     ┊         ┊
    ┊   ┏━┻━┓ ┊   ┃     ┊         ┊
0.20┊   ┃   7 ┊   ┃     ┊   7     ┊
    ┊   ┃   ┃ ┊   ┃     ┊ ┏━┻━┓   ┊
0.17┊   6   ┃ ┊   6     ┊ ┃   ┃   ┊
    ┊ ┏━┻┓  ┃ ┊ ┏━┻━┓   ┊ ┃   ┃   ┊
0.09┊ ┃  5  ┃ ┊ ┃   5   ┊ ┃   5   ┊
    ┊ ┃ ┏┻┓ ┃ ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊
0.07┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊
    ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊
0.00┊ 0 1 3 2 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊
  0.00      2.00      7.00      10.00
*/
char const* unary_ex_nodes = "1  0       0  -1\n"
                             "1  0       0  -1\n"
                             "1  0       0  -1\n"
                             "1  0       0  -1\n"
                             "0  0.071   0  -1\n"
                             "0  0.090   0  -1\n"
                             "0  0.170   0  -1\n"
                             "0  0.202   0  -1\n"
                             "0  0.253   0  -1\n";
char const* unary_ex_edges = "2 10 4 2,3\n"
                             "0 10 5 1\n"
                             "0 2  5 3\n"
                             "2 10 5 4\n"
                             "0 7  6 0,5\n"
                             "7 10 7 0\n"
                             "0 2  7 2\n"
                             "7 10 7 5\n"
                             "0 7  8 6\n"
                             "0 2  8 7\n";

/* We make one mutation for each tree, over unary nodes if they exist */
char const* unary_ex_sites     = "1.0    0\n"
                                 "4.5    0\n"
                                 "8.5    0\n";
char const* unary_ex_mutations = "0    2   1\n"
                                 "1    6   1\n"
                                 "2    5   1\n";

/* An example of a simple tree sequence with multiple marginal trees. */

/* Simple single tree example. */
char const* multiple_tree_ex_nodes = /*                                    */
    "1  0   -1   -1\n"               /*         6        |                */
    "1  0   -1   -1\n"               /*        / \       |                 */
    "1  0   -1   -1\n"               /*       /   \      |     5           */
    "0  1   -1   -1\n"               /*      4     \     |    / \          */
    "0  2   -1   -1\n"               /*     / \     \    |   /   3         */
    "0  3   -1   -1\n"               /*    /   \     \   |  /   / \        */
    "0  4   -1   -1\n";              /*   0     1     2  | 0   1   2       */
                                     /* |----------------|---------------| */
                                     /* 0                1               2 */

char const* multiple_tree_ex_edges = "0.75  1.0   3   1,2\n"
                                     "0.0  0.75   4   0,1\n"
                                     "0.75  1.0   5   0,3\n"
                                     "0.0  0.75   6   2,4\n";

/* Odd topology -- different roots. */

char const* odd_tree1_ex_nodes = /*        |       |   5    */
    "1  0   -1  -1\n"            /*        |   4   |   |    */
    "1  0   -1  -1\n"            /*    3   |   |   |   |    */
    "0  1   -1  -1\n"            /*    |   |   |   |   |    */
    "0  2   -1   -1\n"           /*    2   |   2   |   2    */
    "0  3   -1   -1\n"           /*   / \  |  / \  |  / \   */
    "0  4   -1   -1\n";          /*  0   1 | 0   1 | 0   1  */
                                 /* |------|-------|------| */
                                 /* 0.0    0.2     0.7   1.0*/

char const* odd_tree1_ex_edges = "0.0   1.0 2   0,1\n"
                                 "0.0   0.2 3   2\n"
                                 "0.2   0.7 4   2\n"
                                 "0.7   1.0 4   2\n";

/* An example where some samples descend from other samples, and multiple roots */

char const* multi_root_tree_ex_nodes = "1  0   -1  -1\n" /*  4     5 */
                                       "1  0   -1  -1\n" /*  |     | */
                                       "1  1   -1  -1\n" /*  2     3 */
                                       "1  1   -1  -1\n" /*  |     | */
                                       "0  2   -1  -1\n" /*  0     1 */
                                       "0  2   -1  -1\n";

char const* multi_root_tree_ex_edges = "0   1   2   0\n"
                                       "0   1   3   1\n"
                                       "0   1   4   2\n"
                                       "0   1   5   3\n";

/* Examples of tree sequences where samples have different paths to the same ancestor. */

char const* multi_path_tree_ex_nodes = /*       5        |             */
    "1  0   -1  -1\n"                  /*      / \       |             */
    "1  0   -1  -1\n"                  /*     /   4      |     4       */
    "1  0   -1  -1\n"                  /*    /   / \     |    / \      */
    "0  1   -1  -1\n"                  /*   /   /   \    |   3   \     */
    "0  2   -1  -1\n"                  /*  /   /     \   |  / \   \    */
    "0  3   -1  -1\n";                 /* 0   2       1  | 0   2   1   */
                                       /*----------------|------------ */
                                       /*0.0            0.2         1.0*/

char const* multi_path_tree_ex_edges = "0.2 1.0 3   0\n"
                                       "0.2 1.0 3   2\n"
                                       "0.0 1.0 4   1\n"
                                       "0.0 0.2 4   2\n"
                                       "0.2 1.0 4   3\n"
                                       "0.0 0.2 5   0\n"
                                       "0.0 0.2 5   4\n";

char const* multi_path_tree_ex2_nodes = "1  0   -1  -1\n"
                                        "1  0   -1  -1\n"
                                        "0  1   -1  -1\n"
                                        "0  2   -1  -1\n"
                                        "0  3   -1  -1\n";

char const* multi_path_tree_ex2_edges = "0.6 1.0 2   1\n"
                                        "0.0 1.0 3   0\n"
                                        "0.0 0.6 4   1\n"
                                        "0.6 1.0 4   2\n"
                                        "0.0 1.0 4   3\n";

/* An example of a tree sequence with internally sampled nodes. */

/*
1.20┊         ┊     8   ┊         ┊
    ┊         ┊   ┏━┻━┓ ┊         ┊
1.00┊   7     ┊   ┃   ┃ ┊         ┊
    ┊ ┏━┻━┓   ┊   ┃   ┃ ┊         ┊
0.70┊ ┃   ┃   ┊   ┃   ┃ ┊   6     ┊
    ┊ ┃   ┃   ┊   ┃   ┃ ┊ ┏━┻━┓   ┊
0.50┊ ┃   5   ┊   5   ┃ ┊ ┃   5   ┊
    ┊ ┃ ┏━┻┓  ┊  ┏┻━┓ ┃ ┊ ┃ ┏━┻┓  ┊
0.40┊ ┃ ┃  4  ┊  4  ┃ ┃ ┊ ┃ ┃  4  ┊
    ┊ ┃ ┃ ┏┻┓ ┊ ┏┻┓ ┃ ┃ ┊ ┃ ┃ ┏┻┓ ┊
0.20┊ ┃ ┃ ┃ 3 ┊ ┃ ┃ ┃ 3 ┊ ┃ ┃ ┃ 3 ┊
    ┊ ┃ ┃ ┃   ┊ ┃ ┃ ┃   ┊ ┃ ┃ ┃   ┊
0.10┊ ┃ 1 2   ┊ ┃ 2 1   ┊ ┃ 1 2   ┊
    ┊ ┃       ┊ ┃       ┊ ┃       ┊
0.00┊ 0       ┊ 0       ┊ 0       ┊
  0.00      2.00      8.00      10.00
*/

char const* internal_sample_ex_nodes = "1  0.0   0   -1\n"
                                       "1  0.1   0   -1\n"
                                       "1  0.1   0   -1\n"
                                       "1  0.2   0   -1\n"
                                       "0  0.4   0   -1\n"
                                       "1  0.5   0   -1\n"
                                       "0  0.7   0   -1\n"
                                       "0  1.0   0   -1\n"
                                       "0  1.2   0   -1\n";
char const* internal_sample_ex_edges = "2 8  4 0\n"
                                       "0 10 4 2\n"
                                       "0 2  4 3\n"
                                       "8 10 4 3\n"
                                       "0 10 5 1,4\n"
                                       "8 10 6 0,5\n"
                                       "0 2  7 0,5\n"
                                       "2 8  8 3,5\n";
/* We make one mutation for each tree, some above the internal node */
char const* internal_sample_ex_sites     = "1.0    0\n"
                                           "4.5    0\n"
                                           "8.5    0\n";
char const* internal_sample_ex_mutations = "0    2   1\n"
                                           "1    5   1\n"
                                           "2    5   1\n";

/*** An example of a tree sequence with multiple roots. ***/
/*
0.90┊             ┊         11  ┊             ┊
    ┊             ┊         ┏┻┓ ┊             ┊
0.80┊         10  ┊         ┃ ┃ ┊             ┊
    ┊         ┏┻┓ ┊         ┃ ┃ ┊             ┊
0.40┊     9   ┃ ┃ ┊    9    ┃ ┃ ┊     9       ┊
    ┊   ┏━┻┓  ┃ ┃ ┊  ┏━┻━┓  ┃ ┃ ┊   ┏━┻━━┓    ┊
0.30┊   ┃  ┃  ┃ ┃ ┊  ┃   8  ┃ ┃ ┊   ┃    8    ┊
    ┊   ┃  ┃  ┃ ┃ ┊  ┃  ┏┻┓ ┃ ┃ ┊   ┃   ┏┻┓   ┊
0.20┊   ┃  7  ┃ ┃ ┊  7  ┃ ┃ ┃ ┃ ┊   7   ┃ ┃   ┊
    ┊   ┃ ┏┻┓ ┃ ┃ ┊ ┏┻┓ ┃ ┃ ┃ ┃ ┊ ┏━┻┓  ┃ ┃   ┊
0.10┊   ┃ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┃ ┃ ┊ ┃  6  ┃ ┃   ┊
    ┊   ┃ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┃ ┃ ┊ ┃ ┏┻┓ ┃ ┃   ┊
0.00┊ 5 2 3 4 0 1 ┊ 3 4 1 2 0 5 ┊ 4 0 3 1 2 5 ┊
    0             4             8            10
*/
char const* multiroot_ex_nodes = "1  0.0  0  -1\n"
                                 "1  0.0  0  -1\n"
                                 "1  0.0  0  -1\n"
                                 "1  0.0  0  -1\n"
                                 "1  0.0  0  -1\n"
                                 "1  0.0  0  -1\n"
                                 "0  0.1  0  -1\n"
                                 "0  0.2  0  -1\n"
                                 "0  0.3  0  -1\n"
                                 "0  0.4  0  -1\n"
                                 "0  0.8  0  -1\n"
                                 "0  0.9  0  -1\n";
char const* multiroot_ex_edges = "8  10  6   0,3\n"
                                 "0  8   7   3\n"
                                 "0  10  7   4\n"
                                 "8  10  7   6\n"
                                 "4  10  8   1,2\n"
                                 "0  4   9   2\n"
                                 "0  10  9   7\n"
                                 "4  10  9   8\n"
                                 "0  4   10  0,1\n"
                                 "4  8   11  0,5\n";

/* We make one mutation over each root node */
char const* multiroot_ex_sites     = "1.0    0\n"
                                     "2.0    0\n"
                                     "3.0    0\n"
                                     "5.0    0\n"
                                     "6.0    0\n"
                                     "8.0    0\n"
                                     "9.0    0\n";
char const* multiroot_ex_mutations = "0    10  1\n"
                                     "1    9   1\n"
                                     "2    5   1\n"
                                     "3    11  1\n"
                                     "4    9   1\n"
                                     "5    9   1\n"
                                     "6    5   1\n";

/*** An example of a empty tree sequence. ***/
char const* empty_ex_nodes = "1  0.0  0  -1\n"
                             "1  0.0  0  -1\n"
                             "1  0.0  0  -1\n"
                             "1  0.0  0  -1\n"
                             "1  0.0  0  -1\n"
                             "1  0.0  0  -1\n";
char const* empty_ex_edges = "";

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

/* Returns a tree sequence consisting of a single tree with n samples. This
 * is a full example of the data model, with values included for all fields.
 */
tsk_treeseq_t* caterpillar_tree(tsk_size_t n, tsk_size_t num_sites, tsk_size_t num_mutations) {
    int                    ret;
    tsk_id_t               ret_id;
    tsk_treeseq_t*         ts = reinterpret_cast<tsk_treeseq_t*>(tsk_malloc(sizeof(tsk_treeseq_t)));
    tsk_table_collection_t tables;
    tsk_id_t               j, k, last_node, u;
    int                    state, m;
    double                 position[2];
    tsk_id_t               parents[2]         = {-1, -1};
    char const*            states[]           = {"0", "1"};
    char const*            metadata[]         = {"This", "is", "some", "metadata"};
    int const              num_metadatas      = sizeof(metadata) / sizeof(*metadata);
    char const*            metadata_schema    = "mock metadata schema";
    char const*            ts_metadata        = "This is a caterpillar tree";
    char const*            ts_metadata_schema = "The metadata is an example";
    char const*            prov_timestamp     = "a timestamp, should be ISO8601";
    char const*            prov_record        = "Produced by caterpillar_tree for testing purposes";

    KASSERT(ts != nullptr);
    ret = tsk_table_collection_init(&tables, 0);
    KASSERT(ret == 0);

    KASSERT(num_sites > 0u && num_mutations < n - 1);

    tables.sequence_length = 1.0;

    tsk_table_collection_set_metadata(&tables, ts_metadata, strlen(ts_metadata));
    tsk_table_collection_set_metadata_schema(&tables, ts_metadata_schema, strlen(ts_metadata_schema));
    tsk_reference_sequence_set_metadata_schema(
        &tables.reference_sequence,
        ts_metadata_schema,
        strlen(ts_metadata_schema)
    );
    tsk_reference_sequence_set_metadata(&tables.reference_sequence, ts_metadata, strlen(ts_metadata));
    tsk_reference_sequence_set_data(&tables.reference_sequence, "A", 1);
    tsk_reference_sequence_set_url(&tables.reference_sequence, "B", 1);

    tsk_population_table_set_metadata_schema(&tables.populations, metadata_schema, strlen(metadata_schema));
    tsk_individual_table_set_metadata_schema(&tables.individuals, metadata_schema, strlen(metadata_schema));
    tsk_node_table_set_metadata_schema(&tables.nodes, metadata_schema, strlen(metadata_schema));
    tsk_edge_table_set_metadata_schema(&tables.edges, metadata_schema, strlen(metadata_schema));
    tsk_site_table_set_metadata_schema(&tables.sites, metadata_schema, strlen(metadata_schema));
    tsk_mutation_table_set_metadata_schema(&tables.mutations, metadata_schema, strlen(metadata_schema));
    tsk_migration_table_set_metadata_schema(&tables.migrations, metadata_schema, strlen(metadata_schema));

    for (j = 0; j < (tsk_id_t)n; j++) {
        position[0] = j;
        position[1] = j;
        m           = j % num_metadatas;
        ret_id      = tsk_population_table_add_row(&tables.populations, metadata[m], strlen(metadata[m]));
        KASSERT(ret_id == j);
        ret_id = tsk_individual_table_add_row(
            &tables.individuals,
            0,
            position,
            2,
            parents,
            2,
            metadata[m],
            strlen(metadata[m])
        );
        KASSERT(ret_id == j);
        ret_id = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0, j, j, metadata[m], strlen(metadata[m]));
        KASSERT(ret_id == j);
    }
    last_node = 0;
    for (j = 0; asserting_cast<size_t>(j) < n - 1; j++) {
        m      = j % num_metadatas;
        ret_id = tsk_node_table_add_row(
            &tables.nodes,
            0,
            j + 1,
            j % asserting_cast<tsk_id_t>(n),
            TSK_NULL,
            metadata[m],
            strlen(metadata[m])
        );
        KASSERT(ret_id >= 0);
        u      = ret_id;
        ret_id = tsk_edge_table_add_row(&tables.edges, 0, 1, u, last_node, metadata[m], strlen(metadata[m]));
        KASSERT(ret_id >= 0);
        ret_id = tsk_edge_table_add_row(&tables.edges, 0, 1, u, j + 1, metadata[m], strlen(metadata[m]));
        KASSERT(ret_id >= 0);
        last_node = u;
    }
    for (j = 0; static_cast<size_t>(j) < num_sites; j++) {
        m      = j % num_metadatas;
        ret_id = tsk_site_table_add_row(
            &tables.sites,
            (j + 1) / (double)n,
            states[0],
            strlen(states[0]),
            metadata[m],
            strlen(metadata[m])
        );
        KASSERT(ret_id >= 0);
        u     = 2 * asserting_cast<tsk_id_t>(n) - 3;
        state = 0;
        for (k = 0; k < asserting_cast<tsk_id_t>(num_mutations); k++) {
            m      = k % num_metadatas;
            state  = (state + 1) % 2;
            ret_id = tsk_mutation_table_add_row(
                &tables.mutations,
                j,
                u,
                TSK_NULL,
                tables.nodes.time[u],
                states[state],
                strlen(states[state]),
                metadata[m],
                strlen(metadata[m])
            );
            KASSERT(ret_id >= 0);
            u--;
        }
    }
    ret_id = tsk_provenance_table_add_row(
        &tables.provenances,
        prov_timestamp,
        strlen(prov_timestamp),
        prov_record,
        strlen(prov_record)
    );
    KASSERT(ret_id == 0);

    /* TODO make these consistent with the caterpillar tree topology. */
    for (j = 0; j < asserting_cast<tsk_id_t>(n) - 1; j++) {
        m      = j % num_metadatas;
        ret_id = tsk_migration_table_add_row(
            &tables.migrations,
            0,
            1,
            j,
            j,
            j + 1,
            j + 1.5,
            metadata[m],
            strlen(metadata[m])
        );
        KASSERT(ret_id >= 0);
    }

    ret = tsk_table_collection_sort(&tables, 0, 0);
    KASSERT(ret == 0);
    ret = tsk_table_collection_build_index(&tables, 0);
    KASSERT(ret == 0);
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    KASSERT(ret == 0);
    ret = tsk_treeseq_init(ts, &tables, 0);
    KASSERT(ret == 0);
    tsk_table_collection_free(&tables);
    return ts;
}

void unsort_edges(tsk_edge_table_t* edges, size_t start) {
    size_t      j, k;
    size_t      n    = edges->num_rows - start;
    tsk_edge_t* buff = reinterpret_cast<tsk_edge_t*>(tsk_malloc(n * sizeof(tsk_edge_t)));
    KASSERT(buff != nullptr);

    for (j = 0; j < n; j++) {
        k              = start + j;
        buff[j].left   = edges->left[k];
        buff[j].right  = edges->right[k];
        buff[j].parent = edges->parent[k];
        buff[j].child  = edges->child[k];
    }
    for (j = 0; j < n; j++) {
        k                = start + j;
        edges->left[k]   = buff[n - j - 1].left;
        edges->right[k]  = buff[n - j - 1].right;
        edges->parent[k] = buff[n - j - 1].parent;
        edges->child[k]  = buff[n - j - 1].child;
    }
    free(buff);
}

// TODO Delete the functions we don't need
// static int
// tskit_suite_init(void)
// {
//     int fd = -1;
//     static char* template[] = "/tmp/tsk_c_test_XXXXXX";

//     _tmp_file_name = nullptr;
//     _devnull = nullptr;

//     _tmp_file_name = tsk_malloc(sizeof(template));
//     if (_tmp_file_name == nullptr) {
//         return CUE_NOMEMORY;
//     }
//     strcpy(_tmp_file_name, template);
//     fd = mkstemp(_tmp_file_name);
//     if (fd == -1) {
//         return CUE_SINIT_FAILED;
//     }
//     close(fd);
//     _devnull = fopen("/dev/null", "w");
//     if (_devnull == nullptr) {
//         return CUE_SINIT_FAILED;
//     }
//     return CUE_SUCCESS;
// }

// static int
// tskit_suite_cleanup(void)
// {
//     if (_tmp_file_name != nullptr) {
//         unlink(_tmp_file_name);
//         free(_tmp_file_name);
//     }
//     if (_devnull != nullptr) {
//         fclose(_devnull);
//     }
//     return CUE_SUCCESS;
// }

// static void
// handle_cunit_error(void)
// {
//     fprintf(stderr, "CUnit error occured: %d: %s\n", CU_get_error(), CU_get_error_msg());
//     exit(EXIT_FAILURE);
// }

// int
// test_main(CU_TestInfo *tests, int argc, char **argv)
// {
//     int ret;
//     CU_pTest test;
//     CU_pSuite suite;
//     CU_SuiteInfo suites[] = {
//         {
//             .pName = "tskit",
//             .pInitFunc = tskit_suite_init,
//             .pCleanupFunc = tskit_suite_cleanup,
//             .pTests = tests,
//         },
//         CU_SUITE_INFO_NULL,
//     };
//     if (CUE_SUCCESS != CU_initialize_registry()) {
//         handle_cunit_error();
//     }
//     if (CUE_SUCCESS != CU_register_suites(suites)) {
//         handle_cunit_error();
//     }
//     CU_basic_set_mode(CU_BRM_VERBOSE);

//     if (argc == 1) {
//         CU_basic_run_tests();
//     } else if (argc == 2) {
//         suite = CU_get_suite_by_name("tskit", CU_get_registry());
//         if (suite == nullptr) {
//             printf("Suite not found\n");
//             return EXIT_FAILURE;
//         }
//         test = CU_get_test_by_name(argv[1], suite);
//         if (test == nullptr) {
//             printf("Test '%s' not found\n", argv[1]);
//             return EXIT_FAILURE;
//         }
//         CU_basic_run_test(suite, test);
//     } else {
//         printf("usage: %s <test_name>\n", argv[0]);
//         return EXIT_FAILURE;
//     }

//     ret = EXIT_SUCCESS;
//     if (CU_get_number_of_tests_failed() != 0) {
//         printf("Test failed!\n");
//         ret = EXIT_FAILURE;
//     }
//     CU_cleanup_registry();
//     return ret;
// }
