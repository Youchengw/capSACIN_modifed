# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for capSACIN sidecar binary (macOS arm64)."""

import os
import sys
from pathlib import Path

# Ensure the systemSetup directory is on the path
_spec_dir = Path(SPECPATH)
_system_setup = _spec_dir.parent
sys.path.insert(0, str(_system_setup))

a = Analysis(
    [str(_spec_dir / '__main__.py')],
    pathex=[str(_system_setup)],
    binaries=[],
    datas=[
        # Bundle the capsacin library
        (str(_system_setup / 'capsacin'), 'capsacin'),
    ],
    hiddenimports=[
        'MDAnalysis',
        'MDAnalysis.lib',
        'MDAnalysis.lib.formats',
        'MDAnalysis.lib.formats.cython_util',
        'MDAnalysis.lib._transformations',
        'MDAnalysis.coordinates',
        'MDAnalysis.topology',
        'MDAnalysis.analysis',
        'MDAnalysis.analysis.contacts',
        'numpy',
        'numpy.core',
        'numpy.linalg',
        'numpy.random',
        'pandas',
        'pandas.io',
        'scipy',
        'scipy.spatial',
        'scipy.spatial.transform',
        'scipy.linalg',
        'tqdm',
        'mmtf',
        'gridData',
        'networkx',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'matplotlib',
        'PIL',
        'pytest',
        '_pytest',
        'tkinter',
        'IPython',
        'jupyter',
        'notebook',
        'traitlets',
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=None,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=None)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='capsacin-sidecar',
    debug=False,
    bootloader_ignore_signals=False,
    strip=True,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    target_arch='arm64',
    codesign_identity=None,
    entitlements_file=None,
)
