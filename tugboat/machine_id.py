import subprocess as sp

def machine_id() -> str:
    return sp.check_output("hostname").decode('ascii').strip()
