import asyncio
from multiprocessing import freeze_support

from rnalysis.gui import gui

if __name__ == '__main__':
    freeze_support()
    asyncio.run(gui.run())
