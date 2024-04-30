#!/usr/bin/env Rscript

source("experiments/scripts/common.R")
library(argparse)
library(utils)

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
    "-o", "--output-prefix",
    type = "character",
    required = TRUE,
    help = "Output file name (pdf/png/...)"
)
args <- parser$parse_args()

args <- list()
args$input <- "experiments/measurements/datastructure-stats.csv"
args$output_prefix <- "experiments/plots/datastructure-stats"

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
        collection = factor(collection),
        chromosome = factor(chromosome),
    )

style <- c()
style$height <- 120
style$width <- 120
style$unit <- "mm"
style$text_size <- 8
style$legend.key.size <- unit(3, "mm")
style$point_size <- 1

# --- Number of unique subtrees vs number of edge events in tskit ---
ggplot() +
    geom_point(
        data = data %>% filter(variant == "sfkit_dag", variable == "num_unique_subtrees"),
        aes(x = chromosome, y = value, color = collection), shape = 0) +
    geom_point(
        data = data %>% filter(variant == "tskit", variable == "num_edges"),
        aes(x = chromosome, y = value * 2, color = collection), shape = 15) +
    theme_husky(
        style = style,
        axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    gg_eps()
style$height <- 75
style$width <- 150
ggsave("experiments/plots/unique-subtrees-vs-edges.pdf", width = style$width, height = style$height, units = style$unit)
to_jpeg("experiments/plots/unique-subtrees-vs-edges.pdf")

# --- Cumber of overall subtrees compared to number of unique subtrees ---
inner_join(
    data %>% filter(variant == "all", variable %in% c("num_trees", "num_samples")) %>%
        select(-section, -ds_type, -revision, -machine_id, -iteration) %>%
        pivot_wider(names_from = variable, values_from = value) %>%
        mutate(
            # Not counting the leaves
            overall_subtrees = num_trees * (num_samples - 1)
        ) %>%
        select(collection, chromosome, unit, overall_subtrees),
    data %>% filter(variant == "sfkit_dag", variable == "num_unique_subtrees") %>%
        select(collection, chromosome, unit, variable, value) %>%
        pivot_wider(names_from = variable, values_from = value),
    by = c("collection", "chromosome", "unit")
) %>%
    ggplot(
        aes(
            x = chromosome,
            y = num_unique_subtrees / overall_subtrees,
            color = collection,
            shape = collection
        )
    ) +
    geom_point() +
    theme_husky(
        style = style,
        axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    gg_eps()
style$height <- 75
style$width <- 150
ggsave("experiments/plots/unique-subtrees-vs-overall-subtrees.pdf", width = style$width, height = style$height, units = style$unit)
to_jpeg("experiments/plots/unique-subtrees-vs-overall-subtrees.pdf")
