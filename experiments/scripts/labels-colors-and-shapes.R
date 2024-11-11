#!/usr/bin/env Rscript

collection_labels <- c(
    "tgp" = "TGP",
    "sgdp" = "SGDP",
    "unified" = "Unified",
    "simulated_640k" = "Simulated 640k"
)

# Color brewer
# collection_colors <- c(
#     "tgp" = "#ff7f0e",
#     "sgdp" = "#d62728",
#     "unified" = "#1f77b4",
#     "simulated_640k" = "#2ca02c"
# )

# Dominik

dorn_color <- c(
    "blue" = "#377eb8",
    "orange" = "#ff7f00",
    "red" = "#e41a1c",
    "pink" = "#f781bf",
    "bronze" = "#a65628",
    "green" = "#4daf4a",
    "lilac" = "#984ea3",
    "grey" = "#999999",
    "yellow" = "#dede00"
)

collection_colors <- c(
    "tgp" = unname(dorn_color["blue"]),
    "sgdp" = unname(dorn_color["orange"]),
    "unified" = unname(dorn_color["green"]),
    "simulated_640k" = unname(dorn_color["pink"])
)

collection_shapes <- c(
    "tgp" = 0,
    "sgdp" = 3,
    "unified" = 10,
    "simulated_640k" = 18
)

section_labels <- c(
    "compute_subtree_sizes" = "Compute Subtree Sizes",
    "afs" = "AFS",
    "diversity" = "Diversity",
    "tajimas_d" = "Tajima's D",
    "num_segregating_sites" = "Segregating Sites",
    "divergence" = "Divergence",
    "fst" = "Fst",
    "lca_pairwise" = "LCA Pairwise"
)

variant_labels <- c(
    "tskit" = "tskit",
    "gfkit_subtree" = "gfkit subtree",
    "gfkit_bipart" = "gfkit bipartition"
)

short_variant_labels <- c(
    "tskit" = "tskit",
    "gfkit_subtree" = "subtree",
    "gfkit_bipart" = "bipartition"
)

variant_factor <- factor(
    c("tskit", "gfkit_subtree", "gfkit_bipart"),
    levels = c("tskit", "gfkit_subtree", "gfkit_bipart")
)

variant_ds_stat_labels <- c(
    "tskit" = "tskit edge insert + remove",
    "gfkit_subtree" = "gfkit subtrees",
    "gfkit_bipart" = "gfkit bipartitions"
)

variant_shapes <- c(
    "tskit" = 4,
    "gfkit_subtree" = 5,
    "gfkit_bipart" = 16
)

variant_linetypes <- c(
    "tskit" = "solid",
    "gfkit_subtree" = "dotted",
    "gfkit_bipart" = "dashed"
)

lca_labels = c(
    "lca_10th" = "every 10th sample",
    "lca_9th" = "every 9th sample",
    "lca_8th" = "every 8th sample",
    "lca_7th" = "every 7th sample",
    "lca_6th" = "every 6th sample",
    "lca_5th" = "every 5th sample",
    "lca_4th" = "every 4th sample",
    "lca_3th" = "every 3rd sample",
    "lca_2th" = "every 2nd sample"
)

lca_order = c(
    "lca_2th", "lca_3th", "lca_4th", "lca_5th", "lca_6th",
    "lca_7th", "lca_8th", "lca_9th", "lca_10th"
)
