# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['Interface.py'],
    pathex=[],
    binaries=[
        (r'C:\Users\Joushva Kamble\miniconda3\envs\FirstProject\Lib\site-packages\lightgbm\bin\lib_lightgbm.dll', 'lightgbm\lib')   
    ],
    datas=[
        (r'C:\Users\Joushva Kamble\miniconda3\envs\FirstProject\Lib\site-packages\pmapper\smarts_features.txt', 'pmapper'),
        (r'C:\Users\Joushva Kamble\miniconda3\envs\FirstProject\Lib\site-packages\molfeat\data\skey_parameters.csv', 'molfeat/data'),
        (r'C:\Users\Joushva Kamble\miniconda3\envs\FirstProject\Lib\site-packages\datamol\data\salts_solvents.smi', 'datamol/data'),
        ('Image Files', 'Image Files'),
        ('MLModels', 'MLModels')
    ],
    hiddenimports=[
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
        'h5py._hl.files',  
        'dbm.gnu',
        'dbm.ndbm',
        'rdkit.Chem.Pharm2D.SigFactory',
        'molfeat.calc._map4',
        'darkdetect._mac_detect'
    ],
    hookspath=['Main Files/hooks'],
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