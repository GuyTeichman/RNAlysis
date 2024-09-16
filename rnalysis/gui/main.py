import asyncio
import sys
from multiprocessing import freeze_support
from pathlib import Path

if sys.platform == 'win32':
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

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
