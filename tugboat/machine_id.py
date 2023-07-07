import subprocess as sp

def machine_id() -> str:
    return str(sp.check_output("hostname")).strip()
