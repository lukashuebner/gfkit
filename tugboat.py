#!/usr/bin/env python3
import click
from rich.console import Console

from tugboat.cli.ls import cmd_ls
from tugboat.cli.build import cmd_build
from tugboat.cli.clean import cmd_clean
from tugboat.cli.collect import cmd_collect
from tugboat.cli.plot_benchmarks import cmd_plot_benchmarks
from tugboat.cli.plot_tree_stats import cmd_plot_tree_stats
from tugboat.cli.download import cmd_download
from tugboat.cli.compress import cmd_compress
from tugboat.cli.benchmark import cmd_benchmark

from tugboat.log import Logger
from tugboat.config import Config
from tugboat.deps import _Deps
from tugboat.dataset import Datasets
from tugboat.shell import sh_factory

deps = _Deps()
deps.console = Console()
deps.log = Logger(deps.console)
deps.config = Config()
deps.datasets = Datasets.from_csv(deps.config.DATASETS_CSV, deps.config, deps.log)
deps.sh = sh_factory(deps.log)

@click.group()
def cli():
    pass

# --- ls ---
@cli.command()
def ls():
    """List all datasets"""
    cmd_ls()(deps)

# --- build ---
@cli.command()
@click.argument("config", default="release", type=click.Choice(['debug', 'release', 'relwithdeb', 'all']))
def build(config: str):
    """Configure, compile and test the project."""
    if config == "all":
        for config in ['debug', 'release', 'relwithdeb']:
            cmd_build(config)(deps)
    else:
        cmd_build(config)(deps)

# --- clean ---
@cli.command()
@click.option("--ops-bench/--no-ops-bench", default=True)
@click.option("--conversion-bench/--no-conversion-bench", default=False)
@click.option("--forest/--no-forest", default=False)
@click.option("--trees/--no-trees", default=False)
@click.option("--tsz/--no-tsz", default=False)
@click.option("--collected/--no-collected", default=True)
def clean(ops_bench: bool, conversion_bench: bool, forest: bool, trees: bool, tsz: bool, collected: bool):
    """Delete all generated files (default: ops_bench; optional: conversion_bench, forest, trees, tsz, collected)"""
    cmd_clean(ops_bench, conversion_bench, forest, trees, tsz, collected)(deps)

# --- plot-benchmarks ---
@cli.command(name = "plot-benchmarks")
@click.option("--collect/--no-collect", default=True, help="Run --collect before plotting")
def plot_benchmarks(collect: bool):
    """Plot sfkit vs tskit speed benchmarks"""
    if collect:
        cmd_collect()(deps)
    cmd_plot_benchmarks()(deps)

# --- plot-tree-stats ---
@cli.command(name = "plot-tree-stats")
def plot_tree_stats():
    """Plot sfkit vs tskit speed benchmarks"""
    cmd_plot_tree_stats()(deps)

# --- collect ---
@cli.command()
def collect():
    """Collect all the measurement of different datasets into a single .csv file"""
    cmd_collect()(deps)

# --- download ---
@cli.command()
def download():
    """Download all the datasets"""
    cmd_download()(deps)

# --- benchmark ---
@cli.command()
@click.option("--redo", is_flag=True, help="Redo existing measurements?", default=False)
@click.option(
    "--warmup-iterations",
    help="Number of warmup iterations to run",
    default=deps.config.DEFAULT_WARMUP_ITERATIONS
)
@click.option("--iterations", help="Number of iterations to run", default=deps.config.DEFAULT_ITERATIONS)
def benchmark(redo: bool, warmup_iterations: int, iterations: int):
    """Run all benchmarks"""
    cmd_benchmark(redo, warmup_iterations, iterations)(deps)

# --- compress ---
@cli.command()
@click.option("--parallel", help="Number of conversions to run in parallel", default=1)
def compress(parallel: int):
    cmd_compress(parallel)(deps)

# --- program entry point ---
if __name__ == "__main__":
    cli()
