#!/usr/bin/env Rscript

source("experiments/scripts/common.R")
source("experiments/scripts/labels-colors-and-shapes.R")

input_file <- "experiments/measurements/benchmark-results.csv"
# --- Helper functions ---

# The datasets are sorted by chromosome first and collection second.
# dataset_levels <- expand.grid(
#     c("1kg_chr", "sgdp_chr", "unified_chr"),
#     seq(1, 22)
# ) %>%
#     unite("levels", 1:2, sep = "") %>%
#     pull(levels)

# --- Plotting style ---
style <- c()
style$unit <- "mm"
style$text_size <- 10
style$legend.key.size <- unit(3, "mm")
style$point_size <- 1

# --- Load the measurements from the CSV file ---
# TODO Rename 1kg to TGP everywhere (also in generation scripts)
# TODO Rename sfkit to gfkit everywhere (also in gernation scripts)
# TODO Automate collection of results into a single file
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
    mutate(
        section = factor(section),
        variant = factor(variant, levels = levels(variant_factor)),
        iteration = factor(iteration),
        variable = factor(variable)
    ) %>%
    # TODO Fix revision and machine id checks
    select(-revision, -machine_id) %>%
    mutate(
        collection = factor(collection),
        chromosome = factor(chromosome, levels = seq(1, 22))
    )

# Check that the assumptions made in this plotting script are correct.
# If this fails, then the plotting script needs to be updated.
stopifnot(data %>% filter(variable == "runtime") %>% pull(unit) == "ns")
stopifnot(data %>% filter(variable %in% c("virtmem_delta", "rss_delta", "stack_peak", "heap_peak", "heap_delta")) %>% pull(unit) == "byte")

# --- Analysis of the Runtime of sfkit vs tskit ---
runtime_data <- data %>%
    filter(
        variable == "walltime",
    ) %>%
    group_by(section, variant, collection, chromosome, variable) %>%
    rename(walltime_ns = value) %>%
    summarize(
        walltime_ns_median = median(walltime_ns),
        walltime_ns_min = min(walltime_ns),
        walltime_ns_max = max(walltime_ns),
        walltime_ns_q10 = quantile(walltime_ns, 0.1),
        walltime_ns_q90 = quantile(walltime_ns, 0.9),
        walltime_ns_range = walltime_ns_max - walltime_ns_min,
        walltime_ns_range_q90 = walltime_ns_q90 - walltime_ns_q10,
        walltime_ns_sd = sd(walltime_ns),
        walltime_rel_range = walltime_ns_range / walltime_ns_median,
        walltime_rel_range_q90 = walltime_ns_range_q90 / walltime_ns_median,
        walltime_rel_sd = walltime_ns_sd / walltime_ns_median,
        .groups = "drop"
    )

# Quantify the runtime variation
runtime_data %>%
    select(walltime_rel_range, walltime_rel_range_q90, walltime_rel_sd) %>%
    summary()
    

# TODO Rename sfkit_dag to sfkit_bipartition
# TODO De-duplicate code
reference_runtimes <- runtime_data %>%
        filter(variant == "tskit") %>%
        rename(
            ref_walltime_ns_median = walltime_ns_median,
            ref_walltime_ns_min = walltime_ns_min,
            ref_walltime_ns_max = walltime_ns_max,
            ref_walltime_ns_q10 = walltime_ns_q10,
            ref_walltime_ns_q90 = walltime_ns_q90,
        ) %>%
        select(-variant)

speedup_data_bipart <- inner_join(
        reference_runtimes,
        runtime_data %>% filter(variant == "gfkit_bipart"),
        by = c("section", "collection", "chromosome", "variable")
    ) %>%
    mutate(variant = "gfkit_bipart")

speedup_data_subtree <- inner_join(
        reference_runtimes,
        runtime_data %>% filter(variant == "gfkit_subtree"),
        by = c("section", "collection", "chromosome", "variable")
    ) %>%
    mutate(variant = "gfkit_subtree")

speedup_data <- bind_rows(speedup_data_bipart, speedup_data_subtree) %>%
    mutate(
        speedup_median = ref_walltime_ns_median / walltime_ns_median,
        speedup_min = ref_walltime_ns_median / walltime_ns_min,
        speedup_max = ref_walltime_ns_median / walltime_ns_max,
        speedup_q10 = ref_walltime_ns_median / walltime_ns_q10,
        speedup_q90 = ref_walltime_ns_median / walltime_ns_q90,
        .groups = "drop"
    )

# TODO Factor out common code with the  LCA plots
stats_sections <- c(
    "afs", "diversity", "num_segregating_sites", "tajimas_d", "divergence",
    "fst", "f2", "f3", "f4", "lca_pairwise")

max_shown_speedup <- speedup_data %>% filter(section %in% stats_sections) %>% pull(speedup_median) %>% max()
speedup_y_breaks <- seq(1, max_shown_speedup, 1)
speedup_y_limits <- c(1, max_shown_speedup)

# --- Speedup of sfkit over tskit ---
# TODO Unify and re-use code
range_data <- speedup_data %>%
    filter(
        section %in% stats_sections,
        collection %in% c("tgp", "sgdp", "unified", "simulated_640k")
    ) %>%
    group_by(section, collection, variant) %>%
    summarize(
        speedup_min = min(speedup_median),
        speedup_max = max(speedup_median),
        speedup_median = median(speedup_median),
        .groups = "drop"
    )

point_data <- speedup_data %>%
    filter(
        section %in% stats_sections,
        collection %in% c("tgp", "sgdp", "unified", "simulated_640k")
    ) %>%
    group_by(section, collection, variant) %>%
    summarize(
        speedup_min = min(speedup_median),
        speedup_max = max(speedup_median),
        speedup_median = median(speedup_median),
        .groups = "drop"
    )

style$point_size <- 0.75
style$errorbar_width <- style$point_size
style$errorbar_size <- style$errorbar_width
ggplot() +
    geom_errorbar(
        data = range_data,
        mapping = aes(
            x = section,
            y = speedup_median,
            ymin = speedup_min,
            ymax = speedup_max,
            color = collection,
        ),
        size = style$errorbar_size,
        width = style$errorbar_width,
        position = position_dodge2(width = style$point_size)
    ) +
    geom_point(
        data = point_data,
        mapping = aes(
            x = section,
            y = speedup_median,
            color = collection,
        ),
        size = style$point_size,
        position = position_dodge2(width = style$point_size)
    ) +
    # TODO Add error bars
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.25) +
    theme_husky(
        style = style,
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = style$legend.key.size,
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(-9, 0, 0, 0),
    ) +
    scale_y_continuous(limits = speedup_y_limits, breaks = speedup_y_breaks, minor_breaks = NULL) +
    scale_color_manual(values = collection_colors, labels = collection_labels) +
    scale_shape_manual(values = collection_shapes, labels = collection_labels) +
    scale_x_discrete(
        labels = section_labels,
    ) +
    ylab("speedup over tskit") +
    facet_wrap(~variant, scales = "fixed", ncol = 1, labeller = as_labeller(short_variant_labels)) +
    gg_eps()

style$height <- 120
style$width <- 150
ggsave("experiments/plots/speedup.pdf", width = style$width, height = style$height, units = style$unit)
to_jpeg("experiments/plots/speedup.pdf")

point_data %>%
    filter(variant == "gfkit_subtree", section != "lca_pairwise") %>%
    select(speedup_min, speedup_max, speedup_median) %>%
    summary()

point_data %>%
    filter(variant == "gfkit_subtree", section == "lca_pairwise") %>%
    select(speedup_min, speedup_max, speedup_median) %>%
    summary()

point_data %>%
    filter(variant == "gfkit_bipart") %>%
    select(speedup_min, speedup_max, speedup_median) %>%
    summary()


# --- Speedup of sfkit over tskit for LCA stats ---
stats_sections <- c(
    "lca_10th", "lca_9th", "lca_8th", "lca_7th",
    "lca_6th", "lca_5th", "lca_4th", "lca_3th", "lca_2th"
)

max_shown_speedup <- speedup_data %>% filter(section %in% stats_sections) %>% pull(speedup_median) %>% max()
speedup_y_breaks <- seq(1, max_shown_speedup, 1)
speedup_y_limits <- c(1, max_shown_speedup)

# --- Speedup of sfkit over tskit ---
range_data <- speedup_data %>%
    filter(
        section %in% stats_sections,
        collection %in% c("tgp", "sgdp", "unified", "simulated_640k")
    ) %>%
    group_by(section, collection, variant) %>%
    summarize(
        speedup_min = min(speedup_median),
        speedup_max = max(speedup_median),
        speedup_median = median(speedup_median),
        .groups = "drop"
    )

point_data <- speedup_data %>%
    filter(
        section %in% stats_sections,
        collection %in% c("tgp", "sgdp", "unified", "simulated_640k")
    ) %>%
    group_by(section, collection, variant) %>%
    summarize(
        speedup_min = min(speedup_median),
        speedup_max = max(speedup_median),
        speedup_median = median(speedup_median),
        .groups = "drop"
    )

ggplot() +
    geom_errorbar(
        data = range_data,
        mapping = aes(
            x = factor(section, levels = lca_order),
            y = speedup_median,
            ymin = speedup_min,
            ymax = speedup_max,
            color = collection,
        ),
        size = style$errorbar_size,
        width = style$errorbar_width,
        position = position_dodge2(width = style$point_size)
    ) +
    geom_point(
        data = point_data,
        mapping = aes(
            x = factor(section, levels = lca_order),
            y = speedup_median,
            color = collection,
        ),
        size = style$point_size,
        position = position_dodge2(width = style$point_size)
    ) +
    # TODO Add error bars
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.25) +
    theme_husky(
        style = style,
        # legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = style$legend.key.size,
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(-9, 0, 0, 0),
    ) +
    scale_y_log10() +
    scale_color_manual(values = collection_colors, labels = collection_labels) +
    scale_shape_manual(values = collection_shapes, labels = collection_labels) +
    scale_x_discrete(
        labels = lca_labels,
    ) +
    ylab("speedup over tskit") +
    gg_eps()

style$height <- 75
style$width <- 150

ggsave("experiments/plots/speedup-lca.pdf", width = style$width, height = style$height, units = style$unit)
to_jpeg("experiments/plots/speedup-lca.pdf")

point_data %>%
    filter(variant == "gfkit_subtree", section == "lca_10th") %>%
    select(speedup_min, speedup_max, speedup_median) %>%
    summary()

point_data %>%
    filter(variant == "gfkit_subtree", section == "lca_2th") %>%
    select(speedup_min, speedup_max, speedup_median) %>%
    summary()
