
from tugboat.deps import _Deps
from returns.context import RequiresContext

def cmd_plot_tree_stats() -> RequiresContext[_Deps, None]:
    """Re-create the tree stats plots"""
    def factory(deps: _Deps) -> None:   
        sh = deps.sh

        # Extract statistics from our .trees datasets
        sh(
            f"python3 {deps.config.EXTRACT_TREES_FILES_STATS_PY} " +
            ' '.join(deps.datasets.trees_files()) +
            ' > ' + deps.config.TREES_FILES_STATS_CSV
        )
        sh(
            f"Rscript {deps.config.PLOT_TREES_FILES_STATS_R} "
            f"--input {deps.config.TREES_FILES_STATS_CSV} "
            f"--output-prefix {deps.config.SFKIT_VS_TSKIT_BENCH_PLOT_PREFIX}"
        )

    return RequiresContext(factory)
