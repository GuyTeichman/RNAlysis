# RNAlysis-owned PyInstaller hook for the pure-Python ``graphviz`` package.
#
# RNAlysis depends on ``graphviz`` (NOT ``pygraphviz``); it shells out to the ``dot`` executable to
# render ontology DAGs. To keep the frozen app self-contained we bundle the graphviz program
# executables and their plugins/shared libraries.
#
# Historically ``RNAlysis.spec`` reused pyinstaller-hooks-contrib's ``hook-pygraphviz.py`` for this.
# That broke in pyinstaller-hooks-contrib 2026.7, whose ``hook-pygraphviz.py`` collects the graphviz
# executables by running ``import pygraphviz`` in an isolated child process
# (``_get_executables_info()``). Since RNAlysis does not install ``pygraphviz``, the frozen build
# crashed with ``ModuleNotFoundError: No module named 'pygraphviz'``. This standalone hook resolves
# graphviz purely via ``shutil.which('dot')`` (the approach hooks-contrib 2026.6 and earlier used),
# so it never imports pygraphviz and is immune to that upstream drift. ``requirements_pyinstaller``
# is unpinned, so decoupling from that hook is what keeps the freeze reproducible.

import os
import pathlib
import shutil

from PyInstaller import compat
from PyInstaller.depend import bindepend
from PyInstaller.utils.hooks import logger


def _collect_graphviz_files():
    binaries = []
    datas = []

    # The ``graphviz`` package requires the graphviz programs on PATH. Resolve ``dot`` to locate them.
    dot_binary = shutil.which('dot')
    if not dot_binary:
        logger.warning("hook-graphviz: 'dot' program not found in PATH!")
        return binaries, datas
    logger.info("hook-graphviz: found 'dot' program: %r", dot_binary)
    bin_dir = pathlib.Path(dot_binary).parent

    # Graphviz layout/format programs that may be invoked at runtime. On macOS/Linux several of
    # these are symlinks to a single executable.
    progs = (
        "neato", "dot", "twopi", "circo", "fdp", "nop", "osage", "patchwork", "gc",
        "acyclic", "gvpr", "gvcolor", "ccomps", "sccmap", "tred", "sfdp", "unflatten",
    )

    logger.debug("hook-graphviz: collecting graphviz program executables...")
    for program_name in progs:
        program_binary = shutil.which(program_name)
        if not program_binary:
            logger.debug("hook-graphviz: graphviz program %r not found!", program_name)
            continue
        # Only collect programs from the same directory as ``dot`` so we don't pull in an unrelated
        # graphviz install that happens to be on PATH.
        if pathlib.Path(program_binary).parent != bin_dir:
            logger.debug(
                "hook-graphviz: found program %r (%r) outside of directory %r - ignoring!",
                program_name, program_binary, str(bin_dir))
            continue
        logger.debug("hook-graphviz: collecting graphviz program %r: %r", program_name, program_binary)
        binaries += [(program_binary, '.')]

    # The graphviz shared libraries are picked up automatically by PyInstaller's binary-dependency
    # analysis of the collected executables; we still have to collect the plugins + their config file.
    logger.debug("hook-graphviz: looking for graphviz plugin directory...")
    if compat.is_win:
        # On Windows every install variant keeps the plugins and config next to the executables.
        plugin_dir = bin_dir
        plugin_dest_dir = '.'  # Collect into the top-level application directory.
        plugin_pattern = '*gvplugin*.dll'
    else:
        # Resolve the graphviz shared-library directory from ``dot``'s binary imports, then find the
        # ``graphviz`` plugin sub-directory next to it.
        graphviz_lib_candidates = ['cdt', 'gvc', 'cgraph']
        if hasattr(bindepend, 'get_imports'):
            # PyInstaller >= 6.0
            dot_imports = [path for name, path in bindepend.get_imports(dot_binary) if path is not None]
        else:
            # PyInstaller < 6.0
            dot_imports = bindepend.getImports(dot_binary)

        graphviz_lib_paths = [
            path for path in dot_imports
            if any(candidate in os.path.basename(path) for candidate in graphviz_lib_candidates)
        ]
        if not graphviz_lib_paths:
            logger.warning("hook-graphviz: could not determine location of graphviz shared libraries!")
            return binaries, datas

        graphviz_lib_dir = pathlib.Path(graphviz_lib_paths[0]).parent
        logger.debug("hook-graphviz: location of graphviz shared libraries: %r", str(graphviz_lib_dir))

        plugin_dir = graphviz_lib_dir / 'graphviz'
        plugin_dest_dir = 'graphviz'  # Collect into a graphviz sub-directory.
        if compat.is_darwin:
            plugin_pattern = '*gvplugin*.dylib'
        else:
            # Only the versioned .so plugin files; the unversioned ones are dev-package link targets.
            plugin_pattern = '*gvplugin*.so.*'

    if not plugin_dir.is_dir():
        logger.warning("hook-graphviz: could not determine location of graphviz plugins!")
        return binaries, datas

    logger.info("hook-graphviz: collecting graphviz plugins from directory: %r", str(plugin_dir))
    binaries += [(str(file), plugin_dest_dir) for file in plugin_dir.glob(plugin_pattern)]
    datas += [(str(file), plugin_dest_dir) for file in plugin_dir.glob("config*")]  # e.g. ``config6``

    return binaries, datas


binaries, datas = _collect_graphviz_files()
