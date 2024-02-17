import asyncio
from multiprocessing import freeze_support

import nest_asyncio

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
