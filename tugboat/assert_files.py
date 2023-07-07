from os.path import isfile, getsize, exists
from rich.console import Console
import sys

def assert_isfile(file: str, console: Console, critical: bool = False) -> None:
    if not isfile(file):
        console.print(f"{file} is not a file")
        if critical:
            sys.exit(1)

def assert_exists(file: str, critical = False) -> None:
    if not exists(file):
        console.log(f"{file} already exists")
        if critical:
            sys.exit(1)

def assert_not_exists(file: str, critical = False) -> None:    
    if exists(file):
        console.log(f"{file} does not exist")
        if critical:
            sys.exit(1)

def assert_not_empty(file: str, critical = False) -> None:
    if getsize(file) == 0:
        console.log(f"{file} is empty")
        if critical:
            sys.exit(1)
