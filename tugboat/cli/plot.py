import subprocess as sp
from tugboat.deps import _Deps
from returns.context import RequiresContext

def cmd_plot() -> RequiresContext[_Deps, None]:
    """Re-create all plots"""
    def factory(deps: _Deps) -> None:   
        sh = deps.sh

        # Compare speed of tskit vs sfkit
        sh(
            f"Rscript {deps.config.PLOT_SFKIT_VS_TSKIT_BENCH_R} "
            f"--input {deps.config.SFKIT_VS_TSKIT_BENCH_CSV} "
            f"--output {deps.config.SFKIT_VS_TSKIT_BENCH_PLOT} "
        )

        # Extract statistics from our .trees datasets
        # TODO
        sh(
            f"python3 {deps.config.EXTRACT_TREES_FILES_STATS_PY} " +
            ' '.join(deps.datasets.trees_files())
        )
        sh(
            f"Rscript {deps.config.PLOT_TREES_FILES_STATS_R} "
            f"--input {deps.config.TREES_FILES_STATS_CSV} "
            f"--output-prefix {deps.config.SFKIT_VS_TSKIT_BENCH_PLOT_PREFIX}"
        )

    return RequiresContext(factory)
