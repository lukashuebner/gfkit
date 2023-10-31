import fs.path


class Config:
    """Configuration for the this project."""

    # --- defaults ---
    DEFAULT_ITERATIONS = 10
    DEFAULT_WARMUP_ITERATIONS = 1

    # --- directories ---
    DATA_DIR = "data/"
    MEASUREMENTS_DIR = "experiments/measurements/"
    SCRIPTS_DIR = "experiments/scripts/"
    PLOTS_DIR = "experiments/plots/"

    # --- data files ---
    EMPIRICAL_DATASETS_CSV = DATA_DIR + "/empirical-datasets.csv"
    SCALING_DATASETS_CSV = DATA_DIR + "/scaling-datasets.csv"
    SFKIT_VS_TSKIT_BENCH_CSV = fs.path.combine(
        MEASUREMENTS_DIR, "sfkit-vs-tskit-bench.csv")
    TREES_FILES_STATS_CSV = fs.path.combine(
        MEASUREMENTS_DIR, "trees-files-stats.csv")

    # --- datasets collections ---
    EMPIRICAL_COLLECTIONS = ["1kg", "unified", "sgdp"]
    SCALING_COLLECTIONS = ["scaling"]

    # --- plot files ---
    SFKIT_VS_TSKIT_BENCH_PLOT = fs.path.combine(
        PLOTS_DIR, "sfkit-vs-tskit-bench.pdf")
    SFKIT_VS_TSKIT_BENCH_PLOT_PREFIX = fs.path.combine(
        PLOTS_DIR, "trees-files-stats-")

    # --- executable files ---
    SFKIT_BIN = "build/release/benchmarks/sfkit-bench"
    BENCHMARK_TSKITS_PYTHON_ONLY_FUNCS_PY = fs.path.combine(
        SCRIPTS_DIR, "benchmark-tskits-python-only-funcs.py")
    PLOT_SFKIT_VS_TSKIT_BENCH_R = fs.path.combine(
        SCRIPTS_DIR, "plot-sfkit-vs-tskit-bench.R")
    EXTRACT_TREES_FILES_STATS_PY = fs.path.combine(
        SCRIPTS_DIR, "extract-trees-files-stats.py")
    PLOT_TREES_FILES_STATS_R = fs.path.combine(
        SCRIPTS_DIR, "plot-trees-files-stats.R")
