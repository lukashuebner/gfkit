from os.path import isfile, exists
from os import remove

import rich.progress
from returns.context import RequiresContext

from tugboat.git_rev import git_rev
from tugboat.machine_id import machine_id
from tugboat.deps import _Deps


def cmd_benchmark(redo: bool, warmup_iterations: int, iterations: int, collections) -> RequiresContext[_Deps, None]:
    """Run the benchmarks on all the datasets"""
    def factory(deps: _Deps) -> None:
        datasets = list(deps.datasets.by_collection(collections))
        console = deps.console
        config = deps.config
        log = deps.log
        sh = deps.sh

        progress_bar_columns = (
            rich.progress.SpinnerColumn(),
            rich.progress.TextColumn(
                "[progress.description]{task.description}"),
            rich.progress.BarColumn(),
            rich.progress.MofNCompleteColumn()
        )
        with rich.progress.Progress(*progress_bar_columns, console=console) as progress:
            task = progress.add_task(
                "Running benchmarks...", total=len(datasets))

            for dataset in datasets:
                if not isfile(dataset.trees_file()):
                    log.warn(
                        f"{dataset.trees_file()} does not exists or is not a file")
                elif not isfile(dataset.forest_file()):
                    log.warn(
                        f"{dataset.forest_file()} does not exists or is not a file")
                elif not isfile(dataset.bpforest_file()):
                    log.warn(
                        f"{dataset.bpforest_file()} does not exists or is not a file")
                elif not redo and exists(dataset.ops_bench_file()):
                    log.warn(
                        f"{dataset.ops_bench_file()} already exists but --redo not given")
                elif not redo and exists(dataset.ops_bench_file()):
                    log.warn(
                        f"{dataset.tajimasD_bench_file()} already exists but --redo not given")
                else:
                    ok = True
                    # Run sfkit benchmark
                    ret = sh(f"{config.SFKIT_BIN} benchmark"
                             f" --warmup-iterations={warmup_iterations}"
                             f" --iterations={iterations}"
                             f" --forest-file={dataset.forest_file()}"
                             f" --bp-forest-file={dataset.bpforest_file()}"
                             f" --trees-file={dataset.trees_file()}"
                             f" --revision={git_rev()}"
                             f" --machine={machine_id()}"
                             f" > {dataset.ops_bench_file()}"
                             f" 2> /dev/null"
                    )
                    if ret != 0:
                        remove(dataset.ops_bench_file())
                        ok = False

                    # Run benchmark of those functions that are implemented in Python only in tskit
                    sh(f"python3 {config.BENCHMARK_TSKITS_PYTHON_ONLY_FUNCS_PY}"
                       f" --trees-file={dataset.trees_file()}"
                       f" --iterations={iterations}"
                       f" --revision={git_rev()}"
                       f" --machine={machine_id()}"
                       f" > {dataset.tajimasD_bench_file()}"
                    )
                    if ret != 0:
                        remove(dataset.tajimasD_bench_file())
                        ok = False

                    if ok:
                        log.ok(f"Done benchmarking {dataset.basename()}")
                progress.advance(task)
        log.print_num_errors_and_warnings()

    return RequiresContext(factory)
