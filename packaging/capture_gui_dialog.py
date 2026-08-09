#!/usr/bin/env python
"""Capture PNG screenshots of RNAlysis GUI parameter dialogs, straight from the public API.

The GUI is generated from the API by reflection (see CLAUDE.md: the API is the source of
truth, the GUI reflects it). That means the parameter dialog for any public
``Filter``/``FeatureSet``/``fastq``/``enrichment`` function can be rebuilt and grabbed to a
PNG *headlessly* -- no clicking through the running app. Use this to attach before/after
screenshots to any PR that makes a **visible** GUI change: a new/renamed/removed parameter, a
changed or loosened type annotation, a new ``@readable_name``, edited ``:param:`` help text,
etc. (Behind-the-scenes changes with no visible effect on the dialog do not need a shot.)

Usage
-----
    python packaging/capture_gui_dialog.py filtering.CountFilter.filter_low_reads
    python packaging/capture_gui_dialog.py enrichment.FeatureSet.go_enrichment fastq.trim_adapters
    python packaging/capture_gui_dialog.py filtering.CountFilter.filter_low_reads --out screenshots/

Each TARGET is a dotted path resolved from the ``rnalysis`` package
(``<module>.<Class>.<method>`` or ``<module>.<function>``). One PNG per target is written to
``--out`` (default ``./gui_screenshots``), named ``<dotted.target>.png`` -- a name that is
unique per function, which also keeps parallel captures from clobbering each other.

Fonts / headless notes
----------------------
By default the **native** window platform is used so text renders correctly (proven to work
from a headless agent shell on Windows/macOS with a user session). In a CI runner with *no*
display:
  * Linux  -> run under ``xvfb-run`` so a real framebuffer resolves fonts.
  * Qt's bare ``offscreen`` platform (``--offscreen``) renders text as missing-glyph boxes and
    is only useful for a quick layout/geometry check -- never for documentation.

Exit code is the number of targets that failed to capture (0 == all good).
"""
from __future__ import annotations

import argparse
import importlib
import os
import sys
from pathlib import Path

# Mirror rnalysis/gui/gui.py's main function-application window: these parameters are supplied
# by the app (the object being operated on, the parallel backend, etc.), not chosen in the
# dialog, so a faithful screenshot hides them.
DEFAULT_EXCLUDED_PARAMS = {
    'self', 'backend', 'gui_mode', 'legacy_args', 'function_kwargs', 'fname', 'suppress_warnings',
}


def resolve_target(dotted: str):
    """Resolve e.g. 'filtering.CountFilter.filter_low_reads' to the callable it names.

    Imports the longest importable ``rnalysis.<...>`` module prefix, then walks the remaining
    attributes (class, method / nested function).
    """
    parts = dotted.split('.')
    module = None
    consumed = 0
    for i in range(len(parts), 0, -1):
        try:
            module = importlib.import_module('rnalysis.' + '.'.join(parts[:i]))
            consumed = i
            break
        except ModuleNotFoundError:
            continue
    if module is None:
        raise ValueError(f"could not import any 'rnalysis.*' module from target {dotted!r}")
    obj = module
    for attr in parts[consumed:]:
        obj = getattr(obj, attr)
    if not callable(obj):
        raise TypeError(f"target {dotted!r} resolved to a non-callable {type(obj).__name__}")
    return obj


def capture_dialog(dotted: str, out_dir: Path, excluded: set, help_link, app) -> Path:
    """Build the reflection-generated dialog for one function and save it as a PNG."""
    from rnalysis.gui import gui_windows

    func = resolve_target(dotted)
    func_name = getattr(func, 'readable_name', func.__name__)

    win = gui_windows.FuncExternalWindow(func_name, func, help_link, set(excluded))
    win.init_ui()
    win.show()
    app.processEvents()
    app.processEvents()  # let the layout settle before we measure it
    # Grab the inner form widget (parameter group + start/import/export/close buttons), sized to
    # its full natural height, rather than the window. A *window* can't exceed the physical screen
    # height on the native platform, so a long form would always keep a scrollbar and clip params
    # off the bottom; the inner widget has no such limit, so every parameter renders in one shot.
    form = win.scroll_widget
    form.resize(form.sizeHint())
    app.processEvents()
    app.processEvents()  # a second pass lets the form settle at its full size

    out_path = out_dir / f"{dotted}.png"
    pixmap = form.grab()
    ok = pixmap.save(str(out_path))
    win.close()
    if not ok:
        raise RuntimeError(f"QPixmap.save() returned False for {out_path}")
    return out_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Capture PNG screenshots of RNAlysis GUI parameter dialogs from the public API.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('targets', nargs='+', metavar='TARGET',
                        help="dotted path(s), e.g. filtering.CountFilter.filter_low_reads")
    parser.add_argument('--out', default='gui_screenshots', type=Path,
                        help="output directory (created if missing; default: ./gui_screenshots)")
    parser.add_argument('--exclude', default='', metavar='a,b,c',
                        help="extra parameter names to hide, on top of the defaults")
    parser.add_argument('--show-all', action='store_true',
                        help="do not apply the default parameter exclusions (show self/backend/...)")
    parser.add_argument('--help-link', default=None, metavar='URL',
                        help="add the 'Open documentation' link label the real dialog shows")
    parser.add_argument('--offscreen', action='store_true',
                        help="force Qt's offscreen platform (layout only -- text may be tofu)")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.offscreen:
        os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    # else: leave QT_QPA_PLATFORM untouched so the native platform (with real fonts) is used.

    from PyQt6 import QtWidgets  # imported after the platform decision above

    excluded: set = set() if args.show_all else set(DEFAULT_EXCLUDED_PARAMS)
    excluded |= {p.strip() for p in args.exclude.split(',') if p.strip()}

    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)

    failures = 0
    for target in args.targets:
        try:
            path = capture_dialog(target, out_dir, excluded, args.help_link, app)
            print(f"[ok]   {target} -> {path}")
        except Exception as exc:  # keep going; one bad target shouldn't sink the batch
            failures += 1
            print(f"[FAIL] {target}: {type(exc).__name__}: {exc}", file=sys.stderr)

    if failures:
        print(f"\n{failures} of {len(args.targets)} target(s) failed.", file=sys.stderr)
    return failures


if __name__ == '__main__':
    raise SystemExit(main())
