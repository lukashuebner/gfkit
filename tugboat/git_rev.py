import subprocess as sp

def git_rev() -> str:
    return sp.check_output("git rev-parse --short HEAD", shell=True).decode('ascii').strip()
