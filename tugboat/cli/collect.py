from typing import List
from os.path import isfile, getsize
from itertools import chain
import pandas as pd
from tugboat.deps import _Deps
from returns.context import RequiresContext

def cmd_collect():
    """Collect all the benchmark results into a single file"""
    def factory(deps: _Deps) -> None:
        datasets = deps.datasets
        log = deps.log
        config = deps.config

        ops_bench_files = [ds.ops_bench_file() for ds in datasets.all()]
        conversion_bench_files = [ds.conversion_bench_file() for ds in datasets.all()]

        files: List[str] = []
        for file in chain(ops_bench_files, conversion_bench_files):
            if isfile(file):
                if getsize(file) == 0:
                    log.warn(f"{file} is empty")
                else:
                    files.append(file)
            else:
                pass

        merged = pd.concat([pd.read_csv(file, header=0) for file in files])
        merged.to_csv(config.SFKIT_VS_TSKIT_BENCH_CSV, index=False)

    return RequiresContext(factory)
