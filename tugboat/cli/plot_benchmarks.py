from tugboat.deps import _Deps
from returns.context import RequiresContext

def cmd_plot_benchmarks() -> RequiresContext[_Deps, None]:
    """Re-create the benchmark plots"""
    def factory(deps: _Deps) -> None:   
        sh = deps.sh

        # Compare speed of tskit vs sfkit
        sh(
            f"Rscript {deps.config.PLOT_SFKIT_VS_TSKIT_BENCH_R} "
            f"--input {deps.config.SFKIT_VS_TSKIT_BENCH_CSV} "
            f"--output {deps.config.SFKIT_VS_TSKIT_BENCH_PLOT} "
        )

    return RequiresContext(factory)
