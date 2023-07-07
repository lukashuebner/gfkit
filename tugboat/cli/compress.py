from concurrent.futures import ThreadPoolExecutor
from rich.console import Console
from os.path import isfile, exists
from tugboat.git_rev import git_rev
from tugboat.machine_id import machine_id
from tugboat.deps import _Deps
from returns.context import RequiresContext
from rich.console import Console
from tugboat.dataset import Dataset
from tugboat.log import Logger
import subprocess as sp

# TODO
def _compress_and_check_for_error(sfkit: str, dataset: Dataset, console: Console, log: Logger):
    try:
        sh(
            f"compress --warmup-iterations=0 --iterations=1 --forest-output={dataset.forest_file()}"
            f"--trees-file={dataset.trees_file()} revision={git_rev()}, machine={machine_id()}",
            output=dataset.conversion_bench_file(),
            check=True
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

def cmd_compress(parallel: int, dryrun: bool):
    """Compress all the datasets (.trees) to the .forest format"""
    def factory(deps: _Deps) -> None:
        config = deps.config
        console = deps.console
        datasets = deps.datasets
        log = deps.log

        sfkit_compress = f"{config.SFKIT_BIN} compress"
        with console.status("Converting datasets...") as status:
            with ThreadPoolExecutor(max_workers=parallel) as pool:
                for dataset in datasets.all():
                    if isfile(dataset.forest_file()):
                        log.warn(f"{dataset.forest_file()} already exists")
                    elif exists(dataset.conversion_bench_file()):
                        log.warn(f"{dataset.conversion_bench_file()} already exists")
                    else:
                        pool.submit(_compress_and_check_for_error, sfkit_compress, dataset, console)

    return RequiresContext(factory)
