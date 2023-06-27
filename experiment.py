#!/usr/bin/env python3

import sys
import csv
import click
import sh
import fs
import unittest
import os
import os.path
import rich
from rich import print
from rich.pretty import pprint
from rich.console import Console
from rich.table import Table
from rich.progress import Progress
from rich.prompt import Confirm
from glob import glob
from itertools import filterfalse, chain
import pandas as pd
import urllib.request
from concurrent.futures import ThreadPoolExecutor, wait

BENCHMARK_BIN = "build/release/benchmarks/benchmark"
DEFAULT_ITERATIONS = 10
DEFAULT_WARMUP_ITERATIONS = 1
DATA_DIR = "data/"
MEASUREMENTS_DIR = "experiments/measurements"
SCRIPT_DIR = "experiments/scripts"
PLOT_DIR = "experiments/plots"
DATASETS_CSV = DATA_DIR + "/datasets.csv"
SFKIT_VS_TSKIT_BENCH_CSV = fs.path.combine(MEASUREMENTS_DIR, "/sfkit-vs-tskit-bench.csv")
SFKIT_VS_TSKIT_BENCH_PLOT = fs.path.combine(PLOT_DIR, "/sfkit-vs-tskit-bench.pdf")

def error(msg: str):
    print(f"[red]Error !!![/red] {msg}")
    sys.exit(1)

def warn(msg: str):
    print(f"[yellow]Warning ![/yellow] {msg}")

class Dataset:
    def __init__(self, collection: str, chromosome: str, basename: str, tsz_url: str):
        self._collection = collection
        self._chromosome = chromosome
        self._basename = basename
        self._tsz_url = tsz_url

    def collection(self):
        return self._collection

    def chromosome(self):
        return self._chromosome

    def basename(self):
        return self._basename

    def tsz_url(self):
        return self._tsz_url

    def tsz_file(self):
        return fs.path.combine(DATA_DIR, self.basename() + ".tsz")

    def trees_file(self):
        return fs.path.combine(DATA_DIR, self.basename() + ".trees")

    def forest_file(self):
        return fs.path.combine(DATA_DIR, self.basename() + ".forest")

    def ops_bench_file(self):
        return fs.path.combine(MEASUREMENTS_DIR, self.basename() + ".ops_bench.csv")

    def conversion_bench_file(self):
        return fs.path.combine(MEASUREMENTS_DIR, self.basename() + ".conversion_bench.csv")

class Datasets:
    @classmethod
    def from_csv(cls, filename: str):
        try:
            with open(filename) as file:
                datasets = list()
                reader = csv.reader(file, delimiter=',', quotechar='"')
                row = reader.__next__()
                assert(row[0] == "collection")
                assert(row[1] == "chromosome")
                assert(row[2] == "basename")
                assert(row[3] == "tsz_url")
                for row in reader:
                    datasets.append(
                        Dataset(
                            collection=row[0],
                            chromosome=row[1],
                            basename=row[2],
                            tsz_url=row[3]
                        )
                    )
            return cls(datasets)
        except OSError:
            error(f"Cannot open datasets description: {DATASETS_CSV}")

    def __init__(self, datasets: list[str]):
        self._datasets = datasets

    def all(self):
        return self._datasets

    def by_collection(collection: str):
        return filterfalse(lambda ds: ds.collection == collection, self.datasets)

    # TODO Add existing flag
    # TODO Add the other file types
    def trees_files(self):
        return [ds.trees_file() for ds in self._datasets]

    def ops_bench_files(self):
        return [ds.ops_bench_file() for ds in self._datasets]

    def conversion_bench_files(self):
        return [ds.conversion_bench_file() for ds in self._datasets]

    def tsz_files(self):
        return [ds.tsz_file() for ds in self._datasets]

    def tsz_urls(self):
        return [ds.tsz_url() for ds in self._datasets]

def git_rev():
    return sh.git("rev-parse", "--short", "HEAD").strip()

def machine_id():
    return sh.hostname().strip()

def dry_run_cmd(binary: str, *args, **kwargs):
    cmd = binary + ' ' + ' '.join(args)
    redirect = None
    for key, value in kwargs.items():
        if key == "_out":
            redirect = value
        else:
            cmd = cmd + f" --{key}={value}"

    if redirect is not None:
        cmd = cmd + f" > {redirect}"

    print(cmd)

@click.group()
def cli():
    pass

@cli.command()
def ls():
    """List all datasets"""
    datasets = Datasets.from_csv(DATASETS_CSV)

    table = Table(show_header = True, header_style="bold")
    table.add_column("collection")
    table.add_column("chr", justify="right")
    table.add_column("basename")
    table.add_column("source")

    for ds in datasets.all():
        table.add_row(ds.collection, ds.chromosome, ds.basename, ds.tsz_url)

    Console().print(table)

class DatasetStatusTable:
    def __init__(self):
        self._datasets = list()

        self._table = rich.table.Table()
        self._table = Table(show_header = True, header_style="bold")
        self._table.add_column("collection")
        self._table.add_column("chr", justify="right")
        self._table.add_column("basename")
        self._table.add_column("status")

    def append(self, dataset, status):
        self._table.add_row(dataset.collection, dataset.chromosome, dataset.basename, status)

    def __rich_console__(self, console, options):
        return self._table.__rich_console__(console, options)

@cli.command()
@click.option("--redo", is_flag=True, help="Redo existing measurements?", default=False)
@click.option("--warmup-iterations", help="Number of warmup iterations to run", default=DEFAULT_WARMUP_ITERATIONS)
@click.option("--iterations", help="Number of iterations to run", default=DEFAULT_ITERATIONS)
@click.option("--dry-run", help="Print the commands without running them.", default=False)
def benchmark(redo: bool, warmup_iterations: int, iterations: int, dry_run: bool):

    sfkit_bench = lambda : None
    if dry_run:
        sfkit_bench = lambda *args, **kwargs : dry_run_cmd(BENCHMARK_BIN, *args, **kwargs)
    else:
        if os.path.isfile(BENCHMARK_BIN):
            sfkit_bench = sh.Command(BENCHMARK_BIN)
        else:
            error(f"{BENCHMARK_BIN} does not exist (did you compile?)")

    datasets = Datasets.from_csv(DATASETS_CSV)

    num_errors = 0
    num_warnings = 0
    error_msgs = list()
    console = Console(log_path=False)
    #with console.status("[bold green]Running benchmarks...") as status:
    progress_bar_columns = (
        rich.progress.SpinnerColumn(),
        rich.progress.TextColumn("[progress.description]{task.description}"),
        rich.progress.BarColumn(),
        rich.progress.MofNCompleteColumn()
    )
    with Progress(*progress_bar_columns, console=console) as progress:
        task = progress.add_task("Running benchmarks...", total=len(datasets.all()))
        for ds in datasets.all():#, description = "Running benchmarks..."):
            status = None
            if not os.path.isfile(ds.forest_file()):
                num_warnings += 1
                status = f"[yellow]{ds.forest_file()} does not exists or is not a file"
            elif not redo and os.path.exists(ds.ops_bench_file()):
                num_warnings += 1
                status = f"[yellow]{ds.ops_bench_file()} already exists"
            else:
                try:
                    sfkit_bench(
                            f"--warmup-iterations={warmup_iterations}",
                            f"--forest-input={ds.forest_file()}",
                            f"--trees-file={ds.trees_file()}",
                            iterations=iterations,
                            revision=git_rev(),
                            machine=machine_id(),
                            _out=ds.ops_bench_file()
                    )
                    status = "[bold green]OK"
                except sh.ErrorReturnCode as e:
                    num_errors += 1
                    status = f"[bold red]FAILED with exit code {e.exit_code}"
                    error_msgs.append(
                        f"""
                        [bold]{e.full_cmd}[/bold] failed with exit code [bold]{e.exit_code}[/bold].
                        [lightgrey]{e.stdout}[/lightgrey]
                        [red]{e.stderr}[/red]
                        """
                    )

            assert status is not None
            progress.log(f"\[[bold cyan]{ds.basename}[/bold cyan]] {status}")
            progress.advance(task)

    assert num_errors == len(error_msgs)
    if num_warnings > 0:
        print(f"[bold yellow]There were {num_warnings} warnings")
    if num_errors > 0:
        print(f"[bold red]There were {num_errors} erros")
        for msg in error_msgs:
            print(msg)
    if num_warnings == 0 and num_errors == 0:
        print(f"[bold green]Everything okay :smiley:")

@cli.command()
@click.argument("config", default="release", type=click.Choice(['debug', 'release', 'relwithdeb']))
def compile(config: str):
    console = Console()

    with console.status("Configuring..."):
        sh.cmake(preset=config)
    with console.status("Compiling..."):
        sh.cmake("--build", "--parallel", preset=config)
    with console.status("Testing..."):
        sh.ctest(preset=config)

def rm(path):
    print(path)
    for p in glob(path):
        print(p)
        if os.path.exists(path):
            pass
            # os.remove(path)

@cli.command()
@click.option("--ops-bench/--no-ops-bench", default=True)
@click.option("--conversion-bench/--no-conversion-bench", default=False)
@click.option("--forest/--no-forest", default=False)
@click.option("--trees/--no-trees", default=False)
@click.option("--tsz/--no-tsz", default=False)
@click.option("--collected/--no-collected", default=True)
def clean(ops_bench: bool, conversion_bench: bool, forest: bool, trees: bool, tsz: bool, collected: bool):
    datasets = Datasets.from_csv(DATASETS_CSV)

    files : [str] = []
    if ops_bench:
        files.extend(datasets.ops_bench_files())

    if conversion_bench:
        files.extend(datasets.conversion_bench_files())

    if forest:
        files.extend(datasets.forest_files())

    if trees:
        files.extend(datasets.trees_files())

    if tsz:
        files.extend(datasets.tsz_files())

    if collected:
        files.append(SFKIT_VS_TSKIT_BENCH_CSV)

    # We have to iterate over this list twice
    existing_files = list(filter(lambda file : os.path.isfile(file), files))

    print("I'm going to delete the following files:")
    for file in existing_files:
        print(f"    [bold]-[/bold] {file}")
    is_sure = Confirm.ask("Are you sure?")
    if is_sure:
        for file in existing_files:
            os.remove(file)

def collect_datasets():
    datasets = Datasets.from_csv(DATASETS_CSV)
    ops_bench_files = [ds.ops_bench_file() for ds in datasets.all()]
    conversion_bench_files = [ds.conversion_bench_file() for ds in datasets.all()]

    files: [str] = []
    for file in chain(ops_bench_files, conversion_bench_files):
        if os.path.isfile(file):
            if os.path.getsize(file) == 0:
                warn(f"{file} is empty")
            else:
                files.append(file)
        else:
            pass

    merged = pd.concat([pd.read_csv(file) for file in files])
    merged.to_csv(SFKIT_VS_TSKIT_BENCH_CSV, index=False)

@cli.command()
def collect():
    collect_datasets()

@cli.command()
def plot():
    # Collect all measurement results into a single .csv file
    collect_datasets()

    # Compare speed of tskit vs sfkit
    sh.Rscript(
        fs.path.combine(SCRIPT_DIR, "plot-sfkit-vs-tskit-bench.R"),
        input=SFKIT_VS_TSKIT_BENCH_CSV,
        output=SFKIT_VS_TSKIT_BENCH_PLOT
    )

    # Extract statistics from our .trees datasets
    datasets = Datasets.from_csv(DATASETS_CSV)
    sh.Rscript(
        fs.path.combine(SCRIPT_DIR, "extract-trees-files-stats.py"),
        *datasets.trees_files()
    )

class FileDownloader:
    def __init__(self, files, dl_workers, postprocess_workers, postprocess_text, postprocess_fn, console):
        self._dl_workers = dl_workers
        self._postprocess_workers = postprocess_workers
        self._dl_pool = ThreadPoolExecutor(max_workers=self._dl_workers)
        self._post_pool = ThreadPoolExecutor(max_workers=self._postprocess_workers)
        self._dl_futures = []
        self._post_futures = []
        self._postprocess_test = postprocess_text

        progress_bar_columns = (
            rich.progress.SpinnerColumn(),
            rich.progress.TextColumn("[progress.description]{task.description}"),
            rich.progress.BarColumn(),
            rich.progress.MofNCompleteColumn()
        )
        self._progress = Progress(*progress_bar_columns, console=console)
        self._progress.__enter__()
        self._dl_task = self._progress.add_task("Downloading...", total=len(files))
        self._postprocess_task = self._progress.add_task(f"{postprocess_text}...", total=len(files))
        self._dl_futures = [
            self._dl_pool.submit(self._download, url, dest, postprocess_fn)
            for url, dest in files
        ]

    def _download(self, url, dest, postprocess):
        urllib.request.urlretrieve(url, dest)
        self._progress.log(f"Download of {dst} finished")
        self._progress.advance(self._dl_task)
        self._progress.refresh()
        self._post_futures.append(self._post_pool.submit(_postprocess(postprocess, dest)))

    def _postprocess(self, postprocess, dest):
        postprocess(dest)
        self._progress.log(f"{self._postprocess_text} of {dest} finished")
        self._progress.advance(self._postprocess_task)
        self._progress.refresh()

    def wait(self):
        wait(self._dl_futures)
        wait(self._post_futures)
        self._progress.__exit__()

# TODO Rework the download manager to have access to the datasets and thus the trees file name
def tsunzip(tsz_file, trees_file = None):
    sh.tsunzip(decompress=tsz_file)

@cli.command()
def download():
    console = Console(log_path=False)
    datasets = Datasets.from_csv(DATASETS_CSV)
    downloader = FileDownloader(
        files = [(ds.tsz_url(), ds.tsz_file()) for ds in datasets.all()],
        dl_workers = 2,
        postprocess_workers = 4,
        postprocess_text = "Extracting",
        postprocess_fn = tsunzip,
        console = console
    )
    downloader.wait()

if __name__ == "__main__":
    cli()

