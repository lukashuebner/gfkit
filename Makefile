ITERATIONS = 10
WARMUP_ITERATIONS = 1
PLOT_DIR = experiments/plots
SCRIPT_DIR = experiments/scripts
DATA_DIR = experiments/data
TREES_FILES ?= $(wildcard data/*_chr*.trees)
MEASUREMENT_FILES = $(patsubst data/%.trees, $(DATA_DIR)/%.csv, $(TREES_FILES))
BENCHMARK_BIN = build/release/benchmarks/benchmark

.ONESHELL:

.PHONY: all
all: $(BENCHMARK_BIN) benchmark-all plot-all

.PHONY: benchmark-all
benchmark-all: $(DATA_DIR)/sf-vs-ts-speed.csv

.PHONY: plot-all
plot-all: $(PLOT_DIR)/sf-vs-ts-speed.pdf

# We don't list $(DATA_DIR)/sf-vs-ts-speed.csv as a dependency here, as we do not want the benchmarks to be
# run automatically when we're just trying create the plots.
$(PLOT_DIR)/sf-vs-ts-speed.pdf: $(SCRIPT_DIR)/plot-sf-vs-ts-speed.R $(SCRIPT_DIR)/common.R
	$(SCRIPT_DIR)/plot-sf-vs-ts-speed.R --input "$(DATA_DIR)/sf-vs-ts-speed.csv" --output "$@"

# Rules to collect all benchmark results into a single file.
$(DATA_DIR)/sf-vs-ts-speed.csv: $(MEASUREMENT_FILES)
	cat $^ > "$@"
	sed -i '2,$$ { /^algorithm,variant,dataset,iteration,walltime_ns/d }' "$@"

# Using $(BENCHMARK_BIN) as a dependency here, would lead to re-running the benchmarks on every invocation of Make 
# as $(BENCHMARK_BIN) is a phony target.
$(DATA_DIR)/%.csv: data/%.trees
	$(BENCHMARK_BIN) \
		--warmup-iterations $(WARMUP_ITERATIONS) \
		--iterations $(ITERATIONS) \
		--file "$<" \
		--revision $(file < .git-rev) \
		--machine $(shell hostname) \
		> "$@"

.PHONY: .git-rev
.git-rev:
	git rev-parse --short HEAD > .git-rev
	@status=$$(git status --porcelain);
	test "x$${status}" = x || sed -i "s/$$/-dirty/" .git-rev

.PHONY: $(BENCHMARK_BIN)
$(BENCHMARK_BIN): .git-rev
	cmake --preset release
	cmake --build --preset release --parallel
	ctest --preset release

# Extract statistics from our datasets
$(DATA_DIR)/trees-files-stats.csv: $(TREES_FILES) $(SCRIPT_DIR)/extract-trees-files-stats.py
	$(SCRIPT_DIR)/extract-trees-files-stats.py $(TREES_FILES) > "$@"

TREES_FILES_PLOTS = \
	$(PLOT_DIR)/trees-files-stats-num-trees.pdf \
	$(PLOT_DIR)/trees-files-stats-num-samples.pdf \
	$(PLOT_DIR)/trees-files-stats-num-sites.pdf \
	$(PLOT_DIR)/trees-files-stats-sequence-length.pdf \
	$(PLOT_DIR)/trees-files-stats-num-trees-per-chromosome.pdf \
	$(PLOT_DIR)/trees-files-stats-num-samples-per-chromosome.pdf \
	$(PLOT_DIR)/trees-files-stats-num-sites-per-chromosome.pdf \
	$(PLOT_DIR)/trees-files-stats-sequence-length-per-chromosome.pdf

$(TREES_FILES_PLOTS): $(DATA_DIR)/trees-files-stats.csv $(SCRIPT_DIR)/plot-trees-files-stats.R $(SCRIPT_DIR)/common.R
	$(SCRIPT_DIR)/plot-trees-files-stats.R --input "$<" --output "$@"

# Download and extract datasets
# Extract .trees.tsz files
data/%.trees: data/%.trees.tsz
	tsunzip --decompress $< --stdout > $@

# Simons Genome Diversity Project by Kelleher et al.
data/sgdp_chr%.trees.tsz:
	curl --location --progress-bar --output $@ https://zenodo.org/record/3052359/files/sgdp_chr$*.trees.tsz

# 1kg Project by Kelleher et al.
data/1kg_chr%.trees.tsz:
	curl --location --progress-bar --output $@ https://zenodo.org/record/3051855/files/1kg_chr$*.trees.tsz

# Unified Genealogy of Modern and Ancient Genomes by Wohns et al.
data/unified_chr%.trees.tsz:
	curl --location --progress-bar --output $@ https://zenodo.org/record/5495535/files/hgdp_tgp_sgdp_chr$*_q.dated.trees.tsz

# Simulated genomes by Anderson-Tromce
data/anderson_chr%.trees.tsz:
	curl --location --progress-bar --output $@ https://zenodo.org/record/7702392/files/simulated_chrom_$*.ts.tsz

.PHONY: download-sgdp-data
download-sgdp-data: \
	data/sgdp_chr1.trees \
	data/sgdp_chr2.trees \
	data/sgdp_chr3.trees \
	data/sgdp_chr4.trees \
	data/sgdp_chr5.trees \
	data/sgdp_chr6.trees \
	data/sgdp_chr7.trees \
	data/sgdp_chr8.trees \
	data/sgdp_chr9.trees \
	data/sgdp_chr10.trees \
	data/sgdp_chr11.trees \
	data/sgdp_chr12.trees \
	data/sgdp_chr13.trees \
	data/sgdp_chr14.trees \
	data/sgdp_chr15.trees \
	data/sgdp_chr16.trees \
	data/sgdp_chr17.trees \
	data/sgdp_chr18.trees \
	data/sgdp_chr19.trees \
	data/sgdp_chr20.trees \
	data/sgdp_chr21.trees \
	data/sgdp_chr22.trees

.PHONY: download-1kg-data
download-1kg-data: \
	data/1kg_chr1.trees \
	data/1kg_chr2.trees \
	data/1kg_chr3.trees \
	data/1kg_chr4.trees \
	data/1kg_chr5.trees \
	data/1kg_chr6.trees \
	data/1kg_chr7.trees \
	data/1kg_chr8.trees \
	data/1kg_chr9.trees \
	data/1kg_chr10.trees \
	data/1kg_chr11.trees \
	data/1kg_chr12.trees \
	data/1kg_chr13.trees \
	data/1kg_chr14.trees \
	data/1kg_chr15.trees \
	data/1kg_chr16.trees \
	data/1kg_chr17.trees \
	data/1kg_chr18.trees \
	data/1kg_chr19.trees \
	data/1kg_chr20.trees \
	data/1kg_chr21.trees \
	data/1kg_chr22.trees
	
.PHONY: download-unified-data
download-unified-data: \
	data/unified_chr1.trees \
	data/unified_chr2.trees \
	data/unified_chr3.trees \
	data/unified_chr4.trees \
	data/unified_chr5.trees \
	data/unified_chr6.trees \
	data/unified_chr7.trees \
	data/unified_chr8.trees \
	data/unified_chr9.trees \
	data/unified_chr10.trees \
	data/unified_chr11.trees \
	data/unified_chr12.trees \
	data/unified_chr13.trees \
	data/unified_chr14.trees \
	data/unified_chr15.trees \
	data/unified_chr16.trees \
	data/unified_chr17.trees \
	data/unified_chr18.trees \
	data/unified_chr19.trees \
	data/unified_chr20.trees \
	data/unified_chr21.trees \
	data/unified_chr22.trees

.PHONY: download-anderson-data
download-anderson-data: \
	data/anderson_chr1.trees \
	data/anderson_chr2.trees \
	data/anderson_chr3.trees \
	data/anderson_chr4.trees \
	data/anderson_chr5.trees \
	data/anderson_chr6.trees \
	data/anderson_chr7.trees \
	data/anderson_chr8.trees \
	data/anderson_chr9.trees \
	data/anderson_chr10.trees \
	data/anderson_chr11.trees \
	data/anderson_chr12.trees \
	data/anderson_chr13.trees \
	data/anderson_chr14.trees \
	data/anderson_chr15.trees \
	data/anderson_chr16.trees \
	data/anderson_chr17.trees \
	data/anderson_chr18.trees \
	data/anderson_chr19.trees \
	data/anderson_chr20.trees \
	data/anderson_chr21.trees \
	data/anderson_chr22.trees

.PHONY: download-all-data
download-all-data: download-sgdp-data download-1kg-data download-unified-data
