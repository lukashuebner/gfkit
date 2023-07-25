from returns.context import RequiresContext
from tugboat.deps import _Deps

def cmd_build(config: str) -> RequiresContext[_Deps, None]:
    """Configure, compile, and test the project"""
    def factory(deps: _Deps) -> None:
        sh = deps.sh
        sh(f"cmake --preset {config}", exit_on_failure=True)
        sh(f"cmake --build --preset {config} --parallel", exit_on_failure=True)
        sh(f"ctest --preset {config}", exit_on_failure=True)

    return RequiresContext(factory)
