#!/usr/bin/env Rscript

source("experiments/scripts/common.R")
library(argparse)

# --- Parse command line arguments ---
# parser <- ArgumentParser()
#
# by default ArgumentParser will add an help option
# parser$add_argument(
#     "-i", "--input",
#     type = "character",
#     required = TRUE,
#     help = "Input file name (csv)"
# )
# parser$add_argument(
#     "-o", "--output",
#     type = "character",
#     required = TRUE,
#     help = "Output file name (pdf/png/...)"
# )
# args <- parser$parse_args()

# --- For running in RStudio/VSCode/Emacs ---
args <- list()
# args$input <- "experiments/measurements/sfkit-vs-tskit-bench.csv"
args$input <- "experiments/measurements/scaling/scaling.edge.ops.csv"
args$output <- "experiments/plots/sfkit-vs-tskit-bench-scaling.pdf"

# --- Helper functions ---

# Extract the dataset name from the filename; e.g. "data/1kg_chr1.trees" -> "1kg_chr1"
dataset_from_filename <- function(filename) {
    tools::file_path_sans_ext(basename(filename))
}

# Pretty printing of dataset labels
pretty_print_datasets <- function(filename) {
    filename %>%
        str_replace("1kg_chr", "1KG Chr. ") %>%
        str_replace("sgdp_chr", "SGDP Chr. ") %>%
        str_replace("unified_chr", "Unified Chr. ") %>%
        str_replace("anderson_chr", "Anderson Chr. ")
}

collection_labels <- c(
    "1kg" = "1KG",
    "sgdp" = "SGDP",
    "unified" = "Unified"
    # "anderson" = "Anderson (simulated)"
)

collection_colors <- c(
    "1kg" = "#ff7f0e",
    "sgdp" = "#d62728",
    "unified" = "#1f77b4"
    # "anderson" = "#2ca02c"
)

section_colors <- c(
    "load" = "#ff7f0e",
    "build_sf" = "#d62728",
    "compute_subtree_sizes" = "#1f77b4",
    "afs" = "#2ca02c"
)

section_labels <- c(
    "compute_subtree_sizes" = "Compute Subtree Sizes",
    "afs" = "AFS",
    "diversity" = "Diversity",
    "tajimas_d" = "Tajima's D",
    "num_segregating_sites" = "Segregating Sites",
    "divergence" = "Divergence",
    "fst" = "Fst"
)

# The datasets are sorted by chromosome first and collection second.
dataset_levels <- expand.grid(
    c("1kg_chr", "sgdp_chr", "unified_chr", "anderson_chr"),
    seq(1, 22)
) %>%
    unite("levels", 1:2, sep = "") %>%
    pull(levels)

# --- Plotting style ---
style <- c()
style$height <- 60
style$width <- 85
style$unit <- "mm"
style$text_size <- 8
style$legend.key.size <- unit(3, "mm")
style$point_size <- 1

# --- Load the measurements from the CSV file ---
data <- read_csv(
    args$input,
    col_names = TRUE,
    col_types <- cols(
        section = col_character(),
        variant = col_character(),
        ds_type = col_character(),
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
    mutate(
        section = factor(section),
        variant = factor(variant),
        # TODO Move these preprocessing steps out of the plotting script
        # dataset = factor(dataset_from_filename(dataset), levels = dataset_levels),
        iteration = factor(iteration),
        variable = factor(variable)
    ) %>%
    # separate(
    #     col = dataset,
    #     into = c("collection", "chromosome"),
    #     sep = "_chr",
    #     remove = FALSE,
    # ) %>%
    mutate(
        collection = factor(collection),
        chromosome = factor(chromosome, levels = seq(1, 22))
    )

# Check that the assumptions made in this plotting script are correct.
# If this fails, then the plotting script needs to be updated.
stopifnot(data %>% filter(variable == "runtime") %>% pull(unit) == "ns")
stopifnot(data %>% filter(variable %in% c("virtmem_delta", "rss_delta", "stack_peak", "heap_peak", "heap_delta")) %>% pull(unit) == "byte")
stopifnot(data %>% pull(section) %in% c("load_trees_file", "compute_subtree_sizes", "afs", "diversity", "num_segregating_sites", "tajimas_d", "compress_forest_and_sequence", "save_forest_file", "load_forest_file"))

# --- Analysis of the Runtime of sfkit vs tskit ---
runtime_data <- data %>%
    filter(
        variable == "walltime",
    ) %>%
    group_by(section, variant, ds_type, collection, chromosome, revision, machine_id, variable) %>%
    rename(walltime_ns = value) %>%
    summarize(
        walltime_ns_median = median(walltime_ns),
        walltime_ns_min = min(walltime_ns),
        walltime_ns_max = max(walltime_ns),
        walltime_ns_q10 = quantile(walltime_ns, 0.1),
        walltime_ns_q90 = quantile(walltime_ns, 0.9),
        .groups = "drop"
    )

# TODO Rename sfkit_dag to sfkit_edge
speedup_data <- inner_join(
    # Reference (tskit)
    runtime_data %>%
        filter(variant == "tskit") %>%
        rename(
            ref_walltime_ns_median = walltime_ns_median,
            ref_walltime_ns_min = walltime_ns_min,
            ref_walltime_ns_max = walltime_ns_max,
            ref_walltime_ns_q10 = walltime_ns_q10,
            ref_walltime_ns_q90 = walltime_ns_q90,
        ),
    # sfkit
    runtime_data %>% filter(variant == "sfkit_dag"),
    by = c("section", "ds_type", "collection", "chromosome", "revision", "machine_id", "variable")
) %>%
    mutate(
        speedup_median = ref_walltime_ns_median / walltime_ns_median,
        speedup_min = ref_walltime_ns_median / walltime_ns_min,
        speedup_max = ref_walltime_ns_median / walltime_ns_max,
        speedup_q10 = ref_walltime_ns_median / walltime_ns_q10,
        speedup_q90 = ref_walltime_ns_median / walltime_ns_q90,
        .groups = "drop"
    )

speedup_y_breaks <- seq(0, max(speedup_data$speedup_median * 1.05), 1)
speedup_y_limits <- c(0, max(speedup_data$speedup_median * 1.05))

stats_sections <- c("afs", "diversity", "num_segregating_sites", "tajimas_d", "divergence", "fst", "f2", "f3", "f4")

# --- Speedup of sfkit over tskit ---
speedup_data %>%
    filter(section %in% stats_sections) %>%
    select(-ds_type) %>%
    ggplot(
        aes(
            y = speedup_median,
            ymin = speedup_q10,
            ymax = speedup_q90,
            x = factor(collection, levels = c("5k", "10k", "20k", "40k", "80k", "160k", "320k", "640k")),
            color = factor(section, levels = c("f3", "f4", "f2",  "divergence", "afs", "diversity", "num_segregating_sites", "tajimas_d")),
        ),
    ) +
    # geom_point(size = 0.25, position = position_dodge2(width = 0.25)) +
    geom_point(size = 0.25) +
    geom_line(aes(group = section), size = 0.25) +
    # TODO Add error bars
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.25) +
    theme_husky(
        style = style,
        # legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = style$legend.key.size,
        axis.title.x = element_blank(),
        legend.title = element_blank(),
    ) +
    scale_y_continuous(limits = speedup_y_limits, breaks = speedup_y_breaks, minor_breaks = NULL) +
    # scale_color_shape_manual(
    #     color_values = collection_colors,
    #     labels = collection_labels
    # ) +
    scale_color_dark2() +
    scale_x_discrete(
        labels = section_labels,
        # breaks = section_labels,
    ) +
    ylab("speedup over tskit") +
    gg_eps()

style$height <- 75
style$width <- 150
ggsave("experiments/plots/scaling-sfkit-vs-tskit-speedup-aggregated.pdf", width = style$width, height = style$height, units = style$unit)
to_jpeg("experiments/plots/scaling-sfkit-vs-tskit-speedup-aggregated.pdf")

# --- Runtime of sfkit over tskit ---
runtime_data %>%
    filter(section %in% stats_sections) %>%
    select(-dataset) %>%
    ggplot(
        aes(
            y = ns2ms(walltime_ns_median),
            ymin = ns2ms(walltime_ns_q10),
            ymax = ns2ms(walltime_ns_q90),
            x = chromosome,
            color = variant,
            shape = variant
        ),
    ) +
    geom_points_with_errorbars() +
    facet_grid(rows = vars(section), cols = vars(collection), scales = "fixed") +
    theme_husky(
        style = style,
        # legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = style$legend.key.size,
    ) +
    ylab("runtime [ms]") +
    scale_y_continuous() +
    scale_x_discrete(labels = pretty_print_datasets) +
    # scale_color_shape_manual(
    #     color_values = section_colors,
    #     labels = section_labels
    # ) +
    gg_eps()

style$height <- 150
style$width <- 225
ggsave(args$output, width = style$width, height = style$height, units = style$unit)
