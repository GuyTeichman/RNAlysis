import asyncio
import sys
from multiprocessing import freeze_support

from rnalysis.gui import gui

if __name__ == '__main__':
    freeze_support()

    if sys.platform == 'win32':
        asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

    asyncio.run(gui.run())
