import asyncio
from multiprocessing import freeze_support

from rnalysis.gui import gui
from rnalysis.utils import io

if __name__ == '__main__':
    freeze_support()
    io.run_subprocess(['cutadapt', '--help'])
    asyncio.run(gui.run())
