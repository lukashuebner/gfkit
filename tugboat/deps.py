from rich.console import Console
from tugboat.log import Logger
from tugboat.config import Config
from tugboat.dataset import Datasets
from typing import Callable

class _Deps():
    config: Config
    console: Console
    log: Logger
    datasets: Datasets
    sh: Callable[[str, bool], None]
