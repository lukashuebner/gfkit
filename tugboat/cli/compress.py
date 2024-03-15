from concurrent.futures import ThreadPoolExecutor
from os.path import isfile, exists
from tugboat.git_rev import git_rev
from tugboat.machine_id import machine_id
from tugboat.deps import _Deps
from returns.context import RequiresContext
from tugboat.dataset import EmpiricalDataset, SimulatedDataset
from tugboat.log import Logger
import subprocess as sp
from typing import Callable

# TODO


def _compress_and_check_for_error(sfkit_compress: str, dataset, sh: Callable[[str, bool], int], log: Logger):
    try:
        sh(
            f"{sfkit_compress} "
            f"--trees-file={dataset.trees_file()} "
            f"--forest-file={dataset.forest_file()} "
            f"--bp-forest-file={dataset.bpforest_file()} "
            f"--revision={git_rev()} "
            f"--machine={machine_id()} "
            f" > {dataset.conversion_bench_file()}"
        )
        log.ok(f"Converted {dataset.basename()}")
    except sp.CalledProcessError as e:
        # TODO Abstract away execution of subprocesses with logging
        log.raw(
            f"""
            [bold red]FAILED with exit code {e.exit_code}"
            [bold]{e.full_cmd}[/bold] failed with exit code [bold]{e.exit_code}[/bold].
            [lightgrey]{e.stdout}[/lightgrey]
            [red]{e.stderr}[/red]
            """
        )


def cmd_compress(parallel: int, collections) -> RequiresContext[_Deps, None]:
    """Compress all the datasets (.trees) to the .forest format"""
    def factory(deps: _Deps) -> None:
        sh = deps.sh
        config = deps.config
        console = deps.console
        datasets = deps.datasets.by_collection(collections)
        log = deps.log

        sfkit_compress = f"{config.SFKIT_BIN} compress"
        with console.status("Converting datasets...") as status:
            with ThreadPoolExecutor(max_workers=parallel) as pool:
                for dataset in datasets:
                    if isfile(dataset.forest_file()):
                        log.warn(f"{dataset.forest_file()} already exists")
                    elif exists(dataset.conversion_bench_file()):
                        log.warn(
                            f"{dataset.conversion_bench_file()} already exists")
                    else:
                        pool.submit(_compress_and_check_for_error,
                                    sfkit_compress, dataset, sh, log)

    return RequiresContext(factory)
