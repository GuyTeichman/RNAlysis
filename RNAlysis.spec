# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_submodules

hiddenimports = []
hiddenimports += collect_submodules('sklearn')

block_cipher = None


a = Analysis(
    ['rnalysis_app.py'],
    pathex=[],
    binaries=[],
    datas=[('rnalysis/gui/styles','./rnalysis/gui/styles'),
    ('rnalysis/gui/icons','./rnalysis/gui/icons'),
    ('rnalysis/gui/splash.png','./rnalysis/gui'),
    ('rnalysis/gui/logo_small.png','./rnalysis/gui'),
    ('rnalysis/gui/splash_transparent.png','./rnalysis/gui'),
    ('rnalysis/utils/r_templates','./rnalysis/utils/r_templates'),
    ('rnalysis/favicon.ico','./rnalysis'),
    ('rnalysis/utils/uniprot_dataset_abbreviation_dict.json','./rnalysis/utils')],
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={'pygraphviz':{}},
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

