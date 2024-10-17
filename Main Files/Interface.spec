# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['Interface.py'],
    pathex=[],
    binaries=[
        (r'C:\Users\Joushva Kamble\miniconda3\envs\FirstProject\Lib\site-packages\lightgbm\lib\lib_lightgbm.lib', 'lightgbm\\lib\\')
    ],
    datas=[
        (r'C:\Users\Joushva Kamble\miniconda3\envs\FirstProject\Lib\site-packages\lightgbm\VERSION.txt', 'lightgbm'),
        ('Image Files', 'Image Files'),
        ('MLModels', 'MLModels')
    ],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='EGFR QSAR Modeling Software',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    icon='Image Files/Icon.ico'
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='EGFR QSAR Modeling Software'
)
