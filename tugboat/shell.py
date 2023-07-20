import subprocess as sp
import sys
from tugboat.log import Logger

def sh_factory(log):
    def sh(cmd: str, *args, exit_on_failure: bool = False, **kwargs) -> int:
        """Run a shell command, print the output and exit if it fails"""
        try:
            ret = sp.run(cmd, *args, shell=True, check=True, **kwargs)
            return ret.returncode
        except sp.CalledProcessError as e:
            msg = f"[red][bold]{e.cmd}[/bold] failed with exit code [bold]{e.returncode}[/bold].[/red]"
            if e.stdout is not None:
                msg += f"[lightgrey]{e.stdout}[/lightgrey]"
            if e.stderr is not None:
                msg += f"[red]{e.stderr}[/red]"
            log.raw(msg, level=Logger.Level.ERROR)

            exit_on_failure and sys.exit(e.returncode)

    return sh
