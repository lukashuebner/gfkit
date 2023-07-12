#!/usr/bin/env Rscript

source("experiments/scripts/common.R")
library(argparse)

# --- Parse command line arguments ---
parser <- ArgumentParser()

# by default ArgumentParser will add an help option
parser$add_argument(
    "-i", "--input",
    type = "character",
    required = TRUE,
    help = "Input file name (csv)"
)
parser$add_argument(
    "-o", "--output",
    type = "character",
    required = TRUE,
    help = "Output file name (pdf/png/...)"
)
args <- parser$parse_args()

# --- For running in RStudio/VSCode/Emacs ---
# args <- list()
# args$input <- "experiments/measurements/sfkit-vs-tskit-bench.csv"
# args$output <- "experiments/plots/sfkit-vs-tskit-bench.pdf"

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
    "num_segregating_sits" = "Number of Segregating Sites"
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
        dataset = col_character(),
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
        dataset = factor(dataset_from_filename(dataset), levels = dataset_levels),
        iteration = factor(iteration),
        variable = factor(variable)
    ) %>%
    separate(
        col = dataset,
        into = c("collection", "chromosome"),
        sep = "_chr",
        remove = FALSE,
    ) %>%
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
    group_by(section, variant, dataset, collection, chromosome, revision, machine_id, variable) %>%
    rename(walltime_ns = value) %>%
    summarize(
        walltime_ns_median = median(walltime_ns),
        walltime_ns_min = min(walltime_ns),
        walltime_ns_max = max(walltime_ns),
        walltime_ns_q10 = quantile(walltime_ns, 0.1),
        walltime_ns_q90 = quantile(walltime_ns, 0.9),
        .groups = "drop"
    )

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
    runtime_data %>% filter(variant == "sfkit"),
    by = c("section", "dataset", "collection", "chromosome", "revision", "machine_id", "variable")
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

# --- Speedup of sfkit over tskit ---
speedup_data %>%
    filter(section %in% c("afs", "diversity", "num_segregating_sites", "tajimas_d")) %>%
    select(-dataset) %>%
    ggplot(
        aes(
            y = speedup_median,
            ymin = speedup_q10,
            ymax = speedup_q90,
            x = chromosome,
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
    ylab("speedup over tskit") +
    scale_y_continuous() +
    scale_x_discrete( labels = pretty_print_datasets) +
    # scale_color_shape_manual(
    #     color_values = section_colors,
    #     labels = section_labels
    # ) +
    gg_eps()

# --- Runtime of sfkit over tskit ---
runtime_data %>%
    filter(section %in% c("afs", "diversity", "num_segregating_sites", "tajimas_d")) %>%
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
    scale_x_discrete( labels = pretty_print_datasets) +
    # scale_color_shape_manual(
    #     color_values = section_colors,
    #     labels = section_labels
    # ) +
    gg_eps()

style$height <- 150
style$width <- 225
ggsave(args$output, width = style$width, height = style$height, units = style$unit)
