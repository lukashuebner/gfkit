import subprocess as sp

def git_rev() -> str:
    return str(sp.check_output("git rev-parse --short HEAD", shell=True)).strip()
