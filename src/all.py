import importlib
import types

def reloadall(module):
    importlib.reload(module)
    for child in module:
        if isinstance(child, types.ModuleType):
            reloadall(child)
