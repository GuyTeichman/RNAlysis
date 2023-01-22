# -*- mode: python ; coding: utf-8 -*-
import platform

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

block_cipher = None

a = Analysis(
    ['rnalysis_app.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=get_hook_dirs(),
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

if platform.system() == 'Darwin':
    exe_contents = (pyz, a.scripts, a.binaries, a.zipfiles, a.datas, [],)
    exe_kwargs = dict(runtime_tmpdir=None,)
else:
    splash = Splash('rnalysis/gui/splash.png',
                    binaries=a.binaries,
                    datas=a.datas,
                    text_pos=(180, 510),
                    text_size=12,
                    text_color='black')
    exe_contents = (pyz, splash, a.scripts, [],)
    exe_kwargs = dict(exclude_binaries=True,)

exe = EXE(
    *exe_contents,
    **exe_kwargs,
    name='RNAlysis',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    icon='rnalysis/favicon.ico',
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

if platform.system() != 'Darwin':
    coll = COLLECT(
        a.scripts,
        [],
        a.binaries,
        a.zipfiles,
        a.datas,
        strip=False,
        upx=True,
        upx_exclude=[],
        name='RNAlysis-3.3.0',
    )
