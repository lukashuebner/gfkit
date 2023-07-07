from returns.context import RequiresContext
from tugboat.deps import _Deps

def cmd_compile(config: str) -> RequiresContext[_Deps, None]:
    """Configure, compile, and test the project"""
    def factory(deps: _Deps) -> None:
        sh = deps.sh
        sh(f"cmake --preset {config}")
        sh(f"cmake --build --preset {config} --parallel")
        sh(f"ctest --preset {config}")

    return RequiresContext(factory)
