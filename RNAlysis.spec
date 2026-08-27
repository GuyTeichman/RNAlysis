# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.compat import is_darwin
from PyInstaller.utils.hooks import collect_submodules, collect_data_files

datas = []
binaries = []
hiddenimports = []

tmp_collection = collect_data_files('rnalysis')

for item in tmp_collection:
    if 'videos' in item[0]:
        continue
    datas.append(item)

datas += collect_data_files('pyvis')

hiddenimports += collect_submodules('sklearn')
hiddenimports += collect_submodules('cutadapt')
hiddenimports += collect_submodules('pyarrow')
hiddenimports += collect_submodules('numpy')
hiddenimports += collect_submodules('scipy')
hiddenimports += collect_submodules('polars')
hiddenimports += ['matplotlib.backends.backend_pdf', 'matplotlib.backends.backend_svg',
                  'matplotlib.backends.backend_agg', 'matplotlib.backends.backend_pgf',
                  'matplotlib.backends.backend_ps']
# These packages are loaded lazily at runtime (`lazy_loader.load(...)`, see issue #257), which means
# they are imported by *name* and PyInstaller's static analysis cannot see them. Each one is
# imported by its top-level name only -- the attributes RNAlysis reaches for (e.g.
# `hdbscan.dist_metrics`, `matplotlib_venn.layout`) are all pulled in by the package's own __init__,
# and scikit-learn is already covered by collect_submodules('sklearn') above.
hiddenimports += ['seaborn', 'pandas', 'kmedoids', 'hdbscan', 'matplotlib_venn', 'upsetplot']

# UPX-compressing Qt's core libraries (and the CPython DLL) is a known source of slow
# startup (extra decompression on load), occasional crashes, and AV false positives.
# Patterns are matched against binary filenames via pathlib.PurePath.match() (supports '*').
UPX_EXCLUDE = [
    'Qt6*.dll', 'libQt6*.dylib', 'libQt6*.so*',
    'python3*.dll',
]

# PyInstaller's default bootloader manifest (see PyInstaller.utils.win32.winmanifest) declares
# no DPI awareness, so on a scaled display (125%/150%/...) Windows bitmap-stretches the whole
# bootloader process -- including the Splash() screen, which is rendered by the bootloader
# itself before Qt (and Qt's own DPI-awareness call) ever starts. That stretching is exactly
# what makes the splash look oversized and blurry compared to the in-app QSplashScreen, which
# Qt renders after declaring proper DPI awareness at runtime. Embedding this manifest (based on
# PyInstaller's own default, with dpiAware/dpiAwareness added) fixes it at the source instead of
# guessing a shrink factor for splash.png. Windows-only; PyInstaller ignores/warns on other OSes.
WINDOWS_MANIFEST = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<assembly xmlns="urn:schemas-microsoft-com:asm.v1" manifestVersion="1.0">
  <trustInfo xmlns="urn:schemas-microsoft-com:asm.v3">
    <security>
      <requestedPrivileges>
        <requestedExecutionLevel level="asInvoker" uiAccess="false"></requestedExecutionLevel>
      </requestedPrivileges>
    </security>
  </trustInfo>
  <compatibility xmlns="urn:schemas-microsoft-com:compatibility.v1">
    <application>
      <supportedOS Id="{e2011457-1546-43c5-a5fe-008deee3d3f0}"></supportedOS>
      <supportedOS Id="{35138b9a-5d96-4fbd-8e2d-a2440225f93a}"></supportedOS>
      <supportedOS Id="{4a2f28e3-53b9-4441-ba9c-d69d4a4a6e38}"></supportedOS>
      <supportedOS Id="{1f676c76-80e1-4239-95bb-83d0f6d0da78}"></supportedOS>
      <supportedOS Id="{8e0f7a12-bfb3-4fe8-b9a5-48fd50a15a9a}"></supportedOS>
    </application>
  </compatibility>
  <application xmlns="urn:schemas-microsoft-com:asm.v3">
    <windowsSettings>
      <longPathAware xmlns="http://schemas.microsoft.com/SMI/2016/WindowsSettings">true</longPathAware>
      <dpiAware xmlns="http://schemas.microsoft.com/SMI/2005/WindowsSettings">true/pm</dpiAware>
      <dpiAwareness xmlns="http://schemas.microsoft.com/SMI/2016/WindowsSettings">PerMonitorV2, PerMonitor</dpiAwareness>
    </windowsSettings>
  </application>
  <dependency>
    <dependentAssembly>
      <assemblyIdentity type="win32" name="Microsoft.Windows.Common-Controls" version="6.0.0.0" processorArchitecture="*" publicKeyToken="6595b64144ccf1df" language="*"></assemblyIdentity>
    </dependentAssembly>
  </dependency>
</assembly>
"""

# Use RNAlysis's own committed hook for the `graphviz` package (see hooks/hook-graphviz.py),
# resolved relative to the build CWD (the repo root, where `pyinstaller RNAlysis.spec` is run).
# Previously we copied pyinstaller-hooks-contrib's `hook-pygraphviz.py` here, but 2026.7 made that
# hook `import pygraphviz` (which RNAlysis does not install), breaking the frozen build.
hook_path = 'hooks'

a = Analysis(
    ['rnalysis_app.py'],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[hook_path],
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data)

if is_darwin:
    exe_contents = (pyz, a.scripts, a.binaries, a.zipfiles, a.datas, [],)
    exe_kwargs = dict(runtime_tmpdir=None, icon='rnalysis/favicon.icns')
else:
    splash = Splash('rnalysis/gui/splash.png',
                    binaries=a.binaries,
                    datas=a.datas,
                    text_pos=(175, 510),
                    text_font='Calibri',
                    text_size=16,
                    text_color='black',
                    always_on_top=False)
    exe_contents = (pyz, splash, a.scripts, [],)
    exe_kwargs = dict(exclude_binaries=True, icon='rnalysis/favicon.ico', manifest=WINDOWS_MANIFEST)

exe = EXE(
    *exe_contents,
    **exe_kwargs,
    name='RNAlysis',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=UPX_EXCLUDE,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    hide_console='hide-late',
    # Default '_internal' is a generic name that's easy to mistake for another app's support
    # folder if multiple onedir PyInstaller apps are unpacked side by side; scope it to us.
    contents_directory='RNAlysis_internal',
)

if not is_darwin:
    coll = COLLECT(
        exe,
        splash.binaries,
        a.binaries,
        a.zipfiles,
        a.datas,
        strip=False,
        upx=True,
        upx_exclude=UPX_EXCLUDE,
        name='RNAlysis-4.3.0',
    )
