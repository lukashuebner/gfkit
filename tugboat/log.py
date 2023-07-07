from abc import abstractmethod
from rich.console import Console
import sys
from enum import Enum

class Logger:
    class Level(Enum):
        INFO = 1
        OK = 2
        WARN = 3
        ERROR = 4
        CRITICAL = 5

    def __init__(self, console: Console):
        self._console = console

    def raw(self, msg: str, level: Level) -> None:
        self._console.print(msg)
        if level == Logger.Level.ERROR:
            self._num_errors += 1
        elif level == Logger.Level.WARN:
            self._num_warnings += 1
        elif level == Logger.Level.CRITICAL:
            sys.exit(1)

    def info(self, msg: str) -> None:
        self._console.print(self.format_info(msg))

    def ok(self, msg: str) -> None:
        self._console.print(self.format_ok(msg))

    def warn(self, msg: str) -> None:
        self._console.print(self.format_warn(msg))
        self._num_warnings += 1

    def error(self, msg: str) -> None:
        self._console.print(self.format_error(msg))
        self._num_errors += 1

    def critical(self, msg: str) -> None:
        self._console.print(self.format_error(msg))
        sys.exit(1)

    def format_info(msg: str) -> str:
        return msg

    def format_ok(self, msg: str) -> str:
        return f"[green]OK[/green] {msg}"

    def format_error(self, msg: str) -> str:
        return f"[red]Error !!![/red] {msg}"

    def format_warn(self, msg: str) -> str:
        return f"[yellow]Warning ![/yellow] {msg}"

    def format_critical(self, msg: str) -> str:
        return f"[bold red]Critical !!![/bold red] {msg}"

    def print_num_errors(self) -> None:
        self._console.print(f"[bold red]There were {self._num_errors} errors")

    def print_num_warnings(self) -> None:
        self._console.print(f"[bold yellow]There were {self._num_warnings} warnings")

    def print_num_errors_and_warnings(self) -> None:
        if self._num_errors > 0 or self._num_warnings > 0:
            self.print_num_warnings()
            self.print_num_errors()
        else:
            self.ok(f"Everything okay :smiley:")

    _num_errors = 0
    _num_warnings = 0
