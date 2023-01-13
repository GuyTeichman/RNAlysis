import os


def pytest_configure(config):
    os.environ["QT_DEBUG_PLUGINS"] = "1"
