# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_submodules, collect_all

# datas = [('rnalysis/gui/styles', './rnalysis/gui/styles'),
#          ('rnalysis/gui/icons', './rnalysis/gui/icons'),
#          ('rnalysis/gui/splash.png', './rnalysis/gui'),
#          ('rnalysis/gui/logo_small.png', './rnalysis/gui'),
#          ('rnalysis/gui/splash_transparent.png', './rnalysis/gui'),
#          ('rnalysis/favicon.ico', './rnalysis')]
datas = []
binaries = []
hiddenimports = []

tmp_ret = collect_all('rnalysis')

for item in tmp_ret[0]:
    if 'videos' in item[0] or '__pycache__' in item[0] or '.egg-info' in item[0]:
        continue
    datas.append(item)
binaries += tmp_ret[1]
hiddenimports += tmp_ret[2]

hiddenimports += collect_submodules('sklearn')

block_cipher = None

a = Analysis(
    ['rnalysis_app.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={'pygraphviz': {}},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

splash = Splash('rnalysis/gui/splash.png',
                binaries=a.binaries,
                datas=a.datas,
                text_pos=(180, 510),
                text_size=12,
                text_color='black')


exe = EXE(
    pyz,
    splash,
    a.scripts,
    [],
    exclude_binaries=True,
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
coll = COLLECT(
    exe,
    splash.binaries,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='RNAlysis-3.3.0',
)

