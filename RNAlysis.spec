# -*- mode: python ; coding: utf-8 -*-
from pathlib import Path

from PyInstaller.compat import is_darwin
from PyInstaller.utils.hooks import collect_submodules, collect_data_files
from _pyinstaller_hooks_contrib.hooks.stdhooks import get_hook_dirs

datas = []
binaries = []
hiddenimports = []

tmp_ret = collect_data_files('rnalysis')

for item in tmp_ret:
    if 'videos' in item[0]:
        continue
    datas.append(item)

hiddenimports += collect_submodules('sklearn')
hiddenimports += collect_submodules('cutadapt')
hiddenimports += ['matplotlib.backends.backend_pdf', 'matplotlib.backends.backend_svg',
                  'matplotlib.backends.backend_agg', 'matplotlib.backends.backend_pgf',
                  'matplotlib.backends.backend_ps']

with open(Path(get_hook_dirs()[0]).joinpath('hook-pygraphviz.py')) as infile:
    hook_path = Path('hooks')
    hook_path.mkdir()
    with open(hook_path.joinpath('hook-graphviz.py'), 'w') as outfile:
        outfile.write(infile.read())

block_cipher = None

a = Analysis(
    ['rnalysis_app.py'],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[hook_path],
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

if is_darwin:
    exe_contents = (pyz, a.scripts, a.binaries, a.zipfiles, a.datas, [],)
    exe_kwargs = dict(runtime_tmpdir=None, icon='rnalysis/favicon.icns')
else:
    splash = Splash('rnalysis/gui/splash_transparent.png',
                    binaries=a.binaries,
                    datas=a.datas,
                    text_pos=(175, 510),
                    text_font='Calibri',
                    text_size=16,
                    text_color='black',
                    always_on_top=False)
    exe_contents = (pyz, splash, a.scripts, [],)
    exe_kwargs = dict(exclude_binaries=True, icon='rnalysis/favicon.ico')

exe = EXE(
    *exe_contents,
    **exe_kwargs,
    name='RNAlysis',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
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
        upx_exclude=[],
        name='RNAlysis-3.6.2',
    )
