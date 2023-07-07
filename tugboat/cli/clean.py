from os.path import isfile
from os import remove
from rich.prompt import Confirm
from returns.context import RequiresContext
from tugboat.deps import _Deps
from typing import List

def cmd_clean(ops_bench: bool, conversion_bench: bool, forest: bool, trees: bool, tsz: bool, collected: bool) -> RequiresContext[_Deps, None]:
    def factory(deps: _Deps) -> None:
        datasets = deps.datasets
        console = deps.console

        files: List[str] = []
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
            files.append(deps.config.SFKIT_VS_TSKIT_BENCH_CSV)

        # We have to iterate over this list twice
        existing_files = list(filter(lambda file: isfile(file), files))

        if len(existing_files) == 0:
            console.print("Nothing to delete")
        else:
            console.print("I'm going to delete the following files:")
            for file in existing_files:
                console.print(f"    [bold]-[/bold] {file}")

            is_sure = Confirm.ask("Are you sure?")
            if is_sure:
                for file in existing_files:
                    remove(file)

    return RequiresContext(factory)
