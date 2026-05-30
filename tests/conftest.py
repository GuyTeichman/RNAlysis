import os
import sys

# Force headless offscreen rendering for Qt on macOS CI environments
if sys.platform == 'darwin':
    os.environ["QT_QPA_PLATFORM"] = "offscreen"


def pytest_configure(config):
    os.environ["QT_DEBUG_PLUGINS"] = "1"
