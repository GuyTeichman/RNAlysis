import asyncio
import sys
from multiprocessing import freeze_support

if sys.platform == 'win32':
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

import nest_asyncio

try:
    from rnalysis.gui import gui
except Exception as e:  # if running into related to cache, delete the cache and try again
    try:
        from rnalysis.utils.io import clear_cache
        clear_cache()
    except PermissionError as perm_e:
        raise perm_e from e
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
