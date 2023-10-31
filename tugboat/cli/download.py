from rich.table import Table
from tugboat.deps import _Deps
from returns.context import RequiresContext
from os.path import isfile, getsize
import subprocess as sp

# TODO remove? We don't use this anymore?
class DatasetStatusTable:
    def __init__(self):
        self._datasets = list()

        self._table = Table(show_header = True, header_style="bold")
        self._table.add_column("collection")
        self._table.add_column("chr", justify="right")
        self._table.add_column("basename")
        self._table.add_column("status")

    def append(self, dataset, status):
        self._table.add_row(dataset.collection(), dataset.chromosome(), dataset.basename(), status)

    def __rich_console__(self, console, options):
        return self._table.__rich_console__(console, options)

def cmd_download() -> RequiresContext[_Deps, None]:
    """Download all datasets"""
    def factory(deps: _Deps) -> None:
        datasets = deps.datasets
        config = deps.config
        log = deps.log
        sh = deps.sh

        #with console.status("Downloading and extracting datasets...") as status:
        for ds in datasets.by_collection(config.EMPIRICAL_COLLECTIONS):
            if isfile(ds.trees_file()) and getsize(ds.trees_file()) > 0:
                log.warn(f"{ds.trees_file()} already exists, skipping")
            else:
                if isfile(ds.tsz_file()) and getsize(ds.tsz_file()) > 0:
                    log.warn(f"{ds.tsz_file()} already exists, skipping download")
                else:
                    sh(f"curl --location --output {ds.tsz_file()} {ds.tsz_url()} > /dev/null")
                    log.ok(f"Downloaded {ds.basename()}")

                if not isfile(ds.tsz_file()):
                    log.error(f"{ds.tsz_file()} does not exists or is not a file")
                elif getsize(ds.tsz_file()) == 0:
                    log.error(f"{ds.tsz_file()} is empty")
                else:
                    ret = sh(f"tsunzip --decompress {ds.tsz_file()} --stdout > {ds.trees_file()}")
                    if ret == 0:
                        log.ok(f"Extracted {ds.basename()}")
                    else:
                        log.error(f"Failed to extract {ds.basename()}")

    return RequiresContext(factory)
