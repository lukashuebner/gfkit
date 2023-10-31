from rich.table import Table
from returns.context import RequiresContext
from tugboat.deps import _Deps
from rich.pretty import pprint

def cmd_ls() -> RequiresContext[_Deps, None]:
    """List all datasets"""
    def factory(deps: _Deps) -> None:
        table = Table(show_header = True, header_style="bold")
        table.add_column("collection")
        table.add_column("chr", justify="right")
        table.add_column("basename")
        table.add_column("source")

        pprint(deps.datasets.all())
        for ds in deps.datasets.all():
            table.add_row(ds.collection(), ds.chromosome(), ds.basename(), ds.tsz_url())

        deps.console.print(table)

    return RequiresContext(factory)
