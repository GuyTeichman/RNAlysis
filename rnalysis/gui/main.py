from multiprocessing import freeze_support

from rnalysis.gui import gui


def main():
    freeze_support()
    gui.run()


if __name__ == '__main__':
    freeze_support()
    main()
