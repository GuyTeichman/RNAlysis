import asyncio
from multiprocessing import freeze_support
from pathlib import Path

import appdirs
import nest_asyncio

try:
    from rnalysis.gui import gui
except Exception as e:  # if running into related to cache, delete the cache and try again
    try:
        Path(appdirs.user_cache_dir('RNAlysis')).unlink(missing_ok=True)
    except PermissionError:
        raise e
    from rnalysis.gui import gui


def main():
    freeze_support()
    try:
        asyncio.run(gui.run())
    except RuntimeError:
        nest_asyncio.apply()
        asyncio.run(gui.run())


if __name__ == '__main__':
    freeze_support()
    main()
