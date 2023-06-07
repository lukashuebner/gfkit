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
    "-o", "--output",
    type = "character",
    required = TRUE,
    help = "Output file name (pdf/png/...)"
)
args <- parser$parse_args()

# args <- list()
# args$input <- "experiments/data/trees-files-stats.csv"
# args$output <- "experiments/plots/trees-file-stats"

data <- read_csv(
    args$input,
    col_names = TRUE,
    col_types <- cols(
        filename = col_character(),
        collection = col_character(),
        organism  = col_character(),
        chromosome = col_integer(),
        num_trees = col_double(),
        sequence_length = col_double(),
        num_sample_nodes = col_double(),
        total_size = col_double(),
        num_edges = col_double(),
        num_individuals = col_double(),
        num_migrations = col_double(),
        num_mutations = col_double(),
        num_nodes  = col_double(),
        num_populations = col_double(),
        num_provenances = col_double(),
        num_sites = col_double(), 
    )) %>%
    mutate(
        collection = factor(collection),
        organism = factor(organism),
        chromosome = factor(chromosome),
        total_size_hr = hsize(total_size, standard = "IEC"),
    )

stopifnot(nrow(data) == length(unique(data$filename)))
paste("Number of .trees files:", nrow(data))
paste("Collection:", unique(data$collection))
paste("Organism:", unique(data$organism))
paste("Number of chromosomes:", length(unique(data$chromosome)))

style <- c()
style$height <- 40
style$width <- 120
style$unit <- "mm"
style$text_size <- 8
style$legend.key.size <- unit(3, "mm")
style$point_size <- 1

plot_distribution_of <- function(data, y, ylabel, x, output_postfix) {
    x <- enquo(x)
    y <- enquo(y)
    ggplot(data, aes(y = !!y, x = !!x)) +
        geom_boxplot() +
        theme_husky(style) +
        labs(
            x = element_blank(),
            y = ylabel,
        ) +
        coord_flip() +
        gg_eps()

        output <- paste(args$output, output_postfix, ".pdf", sep = "")
        ggsave(output, width = style$width, height = style$height, units = style$unit)
    }

plot_distribution_of(data, y = num_trees, ylabel = "Number of trees", x = collection, output_postfix = "-num-trees")
plot_distribution_of(data, y = sequence_length, ylabel = "Sequence length", x = collection, output_postfix = "-sequence-length")
plot_distribution_of(data, y = num_sample_nodes, ylabel = "Number of samples", x = collection, output_postfix = "-num-samples")
plot_distribution_of(data, y = num_sites, ylabel = "Number of sites", x = collection, output_postfix = "-num-sites")

plot_for_each_chromosome <- function(data, y, ylabel, output_postfix) {
    y <- enquo(y)
    ggplot(data, aes(x = chromosome, y = !!y, color = collection, shape = collection)) +
        geom_point() +
        theme_husky(style) +
        labs(
            x = element_blank(),
            y = ylabel,
        ) +
        gg_eps()

    output <- paste(args$output, output_postfix, ".pdf", sep = "")
    ggsave(output, width = style$width, height = style$height, units = style$unit)
}

plot_for_each_chromosome(data, y = num_trees, ylabel = "Number of trees", output_postfix = "-num-trees-per-chromosome")
plot_for_each_chromosome(data, y = sequence_length, ylabel = "Sequence length", output_postfix = "-sequence-length-per-chromosome")
plot_for_each_chromosome(data, y = num_sample_nodes, ylabel = "Number of samples", output_postfix = "-num-samples-per-chromosome")
plot_for_each_chromosome(data, y = num_sites, ylabel = "Number of sites", output_postfix = "-num-sites-per-chromosome")
