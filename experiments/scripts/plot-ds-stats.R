#!/usr/bin/env Rscript

source("experiments/scripts/common.R")
source("experiments/scripts/labels-colors-and-shapes.R")

library(utils)
library(scales)
library(xtable)

input_file <- "experiments/measurements/ds-stats.csv"

data <- read_csv(
    input_file,
    col_names = TRUE,
    col_types <- cols(
        section = col_character(),
        variant = col_character(),
        collection = col_character(),
        chromosome = col_integer(),
        revision = col_character(),
        machine_id = col_character(),
        iteration = col_integer(),
        variable = col_character(),
        value = col_double(),
        unit = col_character()
    )
) %>%
    select(-iteration, -machine_id, -section) %>%
    mutate(
        collection = factor(collection),
        chromosome = factor(chromosome),
    )

style <- c()
style$unit <- "mm"
style$text_size <- 10
style$legend.key.size <- unit(3, "mm")
style$point_size <- 1

# --- Number of unique subtrees vs number of edge events in tskit ---
tree_events_data <- bind_rows(
    data %>% filter(variant == "tskit", variable == "num_edges") %>% mutate(value = value * 2) %>% select(-variable),
    data %>% filter(variant == "gfkit_subtree", variable == "num_unique_subtrees") %>% select(-variable),
    data %>% filter(variant == "gfkit_bipart", variable == "num_unique_subtrees") %>% select(-variable),
)
value_y_breaks <- seq(0, max(tree_events_data$value * 1.05), 10000000)
value_y_limits <- c(0, max(tree_events_data$value * 1.05))

value_y_break_factor <- 1000000
ggplot(tree_events_data, aes(
    x = chromosome,
    y = value,
    group = paste(collection, variant),
    color = collection,
    shape = variant,
    linetype = variant
    )) +
    geom_point() +
    geom_line() +
    theme_husky(
        style = style,
        legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(-9, -9, -9, -9),
    ) +
    labs(
        x = "Chromosome",
        color = "Dataset",
        shape = NULL,
        linetype = NULL
    ) +
    scale_y_log10(
        n.breaks = 17,
        minor_breaks = NULL,
        labels = function(x) x / value_y_break_factor,
        name = TeX(paste('$\\times 10^', log10(value_y_break_factor), '$', sep = ''))
    ) +
    scale_shape_manual(values = variant_shapes, labels = variant_ds_stat_labels) +
    scale_linetype_manual(values = variant_linetypes, labels = variant_ds_stat_labels) +
    scale_color_manual(values = collection_colors, labels = collection_labels) +
    gg_eps()

style$height <- 90
style$width <- 150
ggsave("experiments/plots/ds-stats.pdf", width = style$width, height = style$height, units = style$unit)
to_jpeg("experiments/plots/ds-stats.pdf")

# --- Number of unique subtrees / unique bipartitions vs number of overall subtrees ---
options("scipen"=10) 
num_with_sd <- function(x, sd) {
    paste0(
        '\\num{',
        round(x, digits = ceiling(-1 * log10(signif(sd, digit = 1)))),
        '(', signif(sd, digits = 1), ')}'
    )
}

data %>%
    mutate(
        variable = case_when(
            variant == "gfkit_bipart" & variable == "num_unique_subtrees" ~ "num_unique_bipartitions",
            variant == "tskit" & variable == "num_edges" ~ "tskit_num_edges",
            TRUE ~ variable
    )) %>%
    filter(
        variant == "tskit" & variable == "tskit_num_edges" |
        variant == "gfkit_subtree" & variable == "num_unique_subtrees" |
        variant == "gfkit_bipart" & variable == "num_unique_bipartitions" |
        variant == "all" & variable %in% c("num_mutations", "num_samples", "num_trees")
    ) %>%
    select(-variant, -revision, -unit) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    mutate(
        num_subtrees = (2 * num_samples - 1) * num_trees,
        num_bipartitions = tskit_num_edges * num_trees,
        proportion_of_unique_subtrees = num_unique_subtrees / num_subtrees,
        # proportion_of_unique_bipartitions = num_unique_bipartitions / num_bipartitions,
        unique_bipartitions_vs_subtrees = num_unique_bipartitions / num_unique_subtrees
    ) %>%
    select(collection, num_samples, num_unique_subtrees, num_subtrees, proportion_of_unique_subtrees, unique_bipartitions_vs_subtrees) %>%
    group_by(collection) %>%
    summarize(
        num_samples = first(num_samples),
        num_subtrees_min = min(num_subtrees),
        num_subtrees_median = median(num_subtrees),
        num_subtrees_max = max(num_subtrees),
        num_subtrees_mean = mean(num_subtrees),
        num_subtrees_sd = sd(num_subtrees),
        num_unique_subtrees_min = min(num_unique_subtrees),
        num_unique_subtrees_median = median(num_unique_subtrees),
        num_unique_subtrees_max = max(num_unique_subtrees),
        num_unique_subtrees_mean = mean(num_unique_subtrees),
        num_unique_subtrees_sd = sd(num_unique_subtrees),
        proportion_of_unique_subtrees_min = min(proportion_of_unique_subtrees),
        proportion_of_unique_subtrees_median = median(proportion_of_unique_subtrees),
        proportion_of_unique_subtrees_max = max(proportion_of_unique_subtrees),
        proportion_of_unique_subtrees_mean = mean(proportion_of_unique_subtrees),
        proportion_of_unique_subtrees_sd = sd(proportion_of_unique_subtrees),
        # proportion_of_unique_bipartitions_min = min(proportion_of_unique_bipartitions),
        # proportion_of_unique_bipartitions_median = median(proportion_of_unique_bipartitions),
        # proportion_of_unique_bipartitions_max = max(proportion_of_unique_bipartitions),
        # proportion_of_unique_bipartitions_mean = mean(proportion_of_unique_bipartitions),
        # proportion_of_unique_bipartitions_sd = sd(proportion_of_unique_bipartitions),
        unique_bipartitions_vs_subtrees_min = min(unique_bipartitions_vs_subtrees),
        unique_bipartitions_vs_subtrees_median = median(unique_bipartitions_vs_subtrees),
        unique_bipartitions_vs_subtrees_max = max(unique_bipartitions_vs_subtrees),
        unique_bipartitions_vs_subtrees_mean = mean(unique_bipartitions_vs_subtrees),
        unique_bipartitions_vs_subtrees_sd = sd(unique_bipartitions_vs_subtrees),
        .groups = "drop"
    ) %>%
    mutate(
        num_samples_text = paste0('\\num{', num_samples, '}'),
        num_subtrees_text = paste0('\\numrange{', num_subtrees_min, '}{', num_subtrees_max, '}'),
        num_unique_subtrees_text = paste0('\\numrange{', num_unique_subtrees_min, '}{', num_unique_subtrees_max, '}'),
        proportion_of_unique_subtrees_text = num_with_sd(proportion_of_unique_subtrees_mean, proportion_of_unique_subtrees_sd),
        unique_bipartitions_vs_subtrees_text = num_with_sd(unique_bipartitions_vs_subtrees_mean, unique_bipartitions_vs_subtrees_sd)
    ) ->
ds_stats_tbl

print(
    xtable(ds_stats_tbl %>%
        select(
            collection, num_unique_subtrees_text
            # collection, num_samples_text, num_subtrees_text,
            # proportion_of_unique_subtrees_text, unique_bipartitions_vs_subtrees_text
    )),
    type = "latex",
    include.rownames = FALSE,
    sanitize.text.function = identity
)
        
# Dataset | num samples | num subtrees | proportion of subtrees that are unique | proportion of bipartitions that are unique | num_unique_bipartitions / num_unique_subtrees

# ggplot(
#     data = proportions,
#         aes(
#             x = collection,
#             y = proportion,
#             color = label,
#             shape = label
#         )
#     ) +
#     geom_point(position = position_dodge2(width = 1)) +
#     # geom_text(data = text_annotations, aes(
#     #     x = collection,
#     #     y = max,
#     #     label = text,
#     # ),
#     #     vjust = 0,
#     #     hjust = -,
#     #     size = 2
#     # ) +
#     theme_husky(
#         style = style,
#         #axis.text.x = element_text(angle = 30, hjust = 1),
#         legend.position = "bottom",
#         legend.margin = margin(-9, -9, -9, -9),
#     ) +
#     labs(
#         x = NULL,
#         color = NULL,
#         shape = NULL
#     ) +
#     scale_x_discrete(labels = collection_labels) +
#     #scale_y_continuous(breaks = value_y_breaks, minor_breaks = NULL, limits = value_y_limits, labels = label_number()) +
#     #scale_shape_manual(values = variant_shapes, labels = variant_labels) +
#     #scale_linetype_manual(values = variant_linetypes, labels = variant_labels) +
#     #scale_color_manual(values = collection_colors, labels = collection_labels) +
#     gg_eps()

# style$height <- 75
# style$width <- 150
# ggsave("experiments/plots/num-uniq-subtrees.pdf", width = style$width, height = style$height, units = style$unit)
# to_jpeg("experiments/plots/num-uniq-subtrees.pdf")

# --- Storage usage ---
# trees | samples | tskit edges | gfkit subtrees | gfkit bipartitions
# TODO Re-use code
node_id_byte <- 4
site_id_byte <- 4
genomic_position_byte <- 4
genomic_state_byte <- 0.25

data %>%
    select(-revision, - unit) %>%
    filter(variable == "num_edges") %>%
    pivot_wider(names_from = variant, values_from = value) %>%
    group_by(collection) %>%
    summarize(
        # num_samples = first(num_samples),
        # num_trees_sum = sum(num_trees),
        # num_gfkit_subtrees_sum = sum(num_unique_subtrees),
        # num_gfkit_bipartitions_sum = sum(num_unique_bipartitions),
        tskit_edges_sum = sum(tskit),
        gfkit_subtree_edges_sum = sum(gfkit_subtree),
        gfkit_bipart_edges_sum = sum(gfkit_bipart),
        tskit_size_MiB = tskit_edges_sum * (2 * node_id_byte + 2 * genomic_position_byte) / 1024 / 1024,
        gfkit_subtree_size_MiB = gfkit_subtree_edges_sum * 2 * node_id_byte / 1024 / 1024,
        gfkit_bipart_size_MiB = gfkit_bipart_edges_sum * 2 * node_id_byte / 1024 / 1024,
        .groups = "drop"
    ) %>%
    select(-tskit_edges_sum, -gfkit_subtree_edges_sum, -gfkit_bipart_edges_sum)

data %>%
    select(-revision, - unit) %>%
    filter(variable %in% c("num_mutations", "num_sites")) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    group_by(collection) %>%
    summarize(
        ancestral_seq_MiB = sum(num_sites) * genomic_state_byte / 1024 / 1024,
        mutations_MiB = sum(num_mutations) * (genomic_state_byte * 2 + site_id_byte + node_id_byte)  / 1024 / 1024,
        seq_MiB = ancestral_seq_MiB + mutations_MiB,
    )

    mutate(
        variable = case_when(
            variant == "gfkit_bipart" & variable == "num_edges" ~ "",
            variant == "tskit" & variable == "num_edges" ~ "tskit_num_edges",
            TRUE ~ variable
    )) %>%
    filter(
        variant == "tskit" & variable == "tskit_num_edges" |
        variant == "gfkit_subtree" & variable == "num_unique_subtrees" |
        variant == "gfkit_bipart" & variable == "num_unique_bipartitions" |
        variant == "all" & variable %in% c("num_mutations", "num_samples", "num_trees")
    ) %>%
    select(-variant, -revision, -unit) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    group_by(collection) %>%
    summarize(
        num_samples = first(num_samples),
        num_trees_sum = sum(num_trees),
        num_gfkit_subtrees_sum = sum(num_unique_subtrees),
        num_gfkit_bipartitions_sum = sum(num_unique_bipartitions),
        num_tskit_edges_sum = sum(tskit_num_edges),
        tskit_size_byte = num_tskit_edges_sum * (2 * node_id_byte + 2 * genomic_position_byte),
        gfkit_subtree_size_byte = num_gfkit_subtrees_sum * (2 * node_id_byte + 2 * genomic_position_byte),
        .groups = "drop"
    ) %>%
    mutate(
        num_subtrees_text = paste0('\\numrange{', num_subtrees_min, '}{', num_subtrees_max, '}'),
        proportion_of_unique_subtrees_text = num_with_sd(proportion_of_unique_subtrees_mean, proportion_of_unique_subtrees_sd),
        unique_bipartitions_vs_subtrees_text = num_with_sd(unique_bipartitions_vs_subtrees_mean, unique_bipartitions_vs_subtrees_sd)
    ) ->
ds_stats_tbl

#  Average re-use of intermediate results
data %>%
    select(-revision, - unit) %>%
    filter(
        variable %in% c("num_edges", "num_unique_subtrees", "num_trees") &&
            variant %in% c("gfkit_subtree", "gfkit_bipart") ||
        variable == "num_trees" && variant == "all"
    ) %>%
    mutate(
        variable = case_when(
            variant == "gfkit_bipart" & variable == "num_edges" ~ "num_edges_bipart",
            variant == "gfkit_subtree" & variable == "num_edges" ~ "num_edges_subtree",
            variant == "gfkit_bipart" & variable == "num_unique_subtrees" ~ "num_unique_biparts",
            variant == "gfkit_subtree" & variable == "num_unique_subtrees" ~ "num_unique_subtrees",
            TRUE ~ variable
    )) %>%
    select(-variant) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    select(-num_samples, -num_sites, -num_mutations, -num_edges) %>%
    group_by(collection) %>%
    summarize(
        num_edges_subtree_sum = sum(num_edges_subtree),
        num_edges_bipart_sum = sum(num_edges_bipart),
        num_unique_subtrees_sum = sum(num_unique_subtrees),
        num_unique_biparts_sum = sum(num_unique_biparts),
        num_trees_sum = sum(num_trees),
        indegree_bipart_mean = mean(num_edges_bipart_sum / (num_unique_biparts_sum - num_trees_sum)),
        indegree_subtree_mean = mean(num_edges_subtree_sum / (num_unique_subtrees_sum - num_trees_sum)),
        .groups = "drop"
    ) %>%
    select(collection, indegree_bipart_mean, indegree_subtree_mean)
