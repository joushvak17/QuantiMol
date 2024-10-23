# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['Interface.py'],
    pathex=[],
    binaries=[],
    datas=[
        (r'C:\Users\Joushva Kamble\miniconda3\envs\FirstProject\Lib\site-packages\pmapper\smarts_features.txt', 'pmapper'),
        ('Image Files', 'Image Files'),
        ('MLModels', 'MLModels')
    ],
    hiddenimports=[
        'molfeat.trans.fp',  # Add any other hidden imports here
        'rdkit',
        'rdkit.Chem',
        'rdkit.Chem.Draw',
        'PIL',
        'joblib',
        'matplotlib.backends.backend_tkagg',
        'matplotlib.figure',
        'ctypes',
        'tkinter',
        'tkinter.filedialog',
        'tkinter.messagebox',
        'tkinter.ttk',
        'numpy',
        'pandas',
        'scipy',
        'sklearn',
        'lightgbm',
        'datamol',
        'pmapper',
        'ctk',
        'h5py._hl.files',  # Add h5py hidden imports
        'dbm.gnu',
        'dbm.ndbm',
        'rdkit.Chem.Pharm2D.SigFactory',
        'molfeat.calc._map4',
        'darkdetect._mac_detect'
    ],
    hookspath=['Main Files/Hooks'],
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
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='Image Files/Icon.ico'
)