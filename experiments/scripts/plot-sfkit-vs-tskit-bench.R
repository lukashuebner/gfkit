#!/usr/bin/env Rscript

source("experiments/scripts/common.R")
library(argparse)

# Command Line Parsing
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

# args <- list()
# args$input <- "experiments/data/sf-vs-ts-speed.csv"
# args$output <- "experiments/plots/sf-vs-ts-speed.pdf"

dataset_from_filename <- function(filename) {
    tools::file_path_sans_ext(basename(filename))
}

style <- c()
style$height <- 60
style$width <- 85
style$unit <- "mm"
style$text_size <- 8
style$legend.key.size <- unit(3, "mm")
style$point_size <- 1

### Data loading ###
data <- read_csv(
    args$input,
    col_names = TRUE,
    col_types <- cols(
        section = col_character(),
        variant = col_character(),
        dataset = col_character(),
        iteration = col_double(),
        variable = col_character()
        value = col_double()
        unit = col_character()
    )
) %>%
    filter(
        section %in% c("compute_subtree_sizes", "compute_afs")
    ) %>%
    mutate(
        section = factor(section),
        variant = factor(variant),
        dataset = factor(dataset_from_filename(dataset)),
        iteration = factor(iteration),
        variable = factor(variable)
    )

dataset_levels <- expand.grid(
        c("tgp_chr", "1kg_chr", "sgdp_chr", "unified_chr", "anderson_chr"),
        seq(1, 22)
    ) %>%
    unite("levels", 1:2, sep = "") %>%
    pull(levels)

# dataset_levels <- levels(data$dataset) %>%
#     str_replace("tgp_chr", "") %>%
#     as.numeric() %>%
#     sort() %>%
#     paste("tgp_chr", ., sep = "")
data$dataset <- factor(data$dataset, levels = dataset_levels)

tskit_data <- data %>%
    filter(
        variant == "tskit",
    ) %>%
    group_by(section, dataset) %>%
    summarize(
        walltime_ms_median = median(walltime_ms),
        .groups = "drop"
    )

sf_data <- right_join(
    tskit_data %>% rename(reference_walltime_ms = walltime_ms_median),
    data %>%
        filter(
            variant == "sf",
            section == "compute_afs"
        ),
    by = c("section", "dataset")
) %>%
    group_by(section, variant, dataset) %>%
    summarize(
        walltime_ms_median = median(walltime_ms),
        walltime_ms_min = min(walltime_ms),
        walltime_ms_max = max(walltime_ms),
        walltime_ms_q10 = quantile(walltime_ms, 0.1),
        walltime_ms_q90 = quantile(walltime_ms, 0.9),
        speedup_median = reference_walltime_ms / walltime_ms_median,
        speedup_min = reference_walltime_ms / walltime_ms_min,
        speedup_max = reference_walltime_ms / walltime_ms_max,
        speedup_q10 = reference_walltime_ms / walltime_ms_q10,
        speedup_q90 = reference_walltime_ms / walltime_ms_q90,
        .groups = "drop"
    )

# Dataset labels
pretty_print_tgp <- function(tgp_filename) {
    tgp_filename %>%
        str_replace("tgp_chr", "TGP Chr. ")
}

section_colors <- c(
    "load" = "#ff7f0e",
    "build_sf" = "#d62728",
    "compute_subtree_sizes" = "#1f77b4",
    "compute_afs" = "#2ca02c"
)

section_labels <- c(
    "load" = "Load SF from TS file",
    "build_sf" = "Build SF from TS",
    "compute_subtree_sizes" = "Compute subtree sizes",
    "compute_afs" = "Compute AFS"
)

x_breaks <- seq(0, max(sf_data$speedup_median * 1.05), 1)
x_limits <- c(0, max(sf_data$speedup_median * 1.05))

ggplot(sf_data) +
    geom_point(
        aes(
            y = speedup_median,
            x = dataset,
            #color = section,
        ),
        position = position_dodge(width = 0.3),
        size = style$point_size
    ) +
    geom_errorbar(
        aes(
            ymin = speedup_q10,
            ymax = speedup_q90,
            x = dataset,
        ),
        width = 0.1
    ) +
    theme_husky(
        style = style,
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = style$legend.key.size,
    ) +
    ylab("speedup over tskit") +
    scale_y_continuous(
        breaks = x_breaks,
        limits = x_limits
    ) +
    scale_x_discrete(
        labels = pretty_print_tgp,
    ) +
    scale_color_shape_manual(
        color_values = section_colors,
        labels = section_labels
    ) +
    gg_eps()

ggsave(args$output, width = style$width, height = style$height, units = style$unit)