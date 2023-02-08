from multiprocessing import freeze_support
import asyncio
from rnalysis.gui import gui


def main():
    freeze_support()
    asyncio.run(gui.run())


if __name__ == '__main__':
    freeze_support()
    main()
