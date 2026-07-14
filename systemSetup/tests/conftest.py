"""capSACIN test fixtures."""
import os
import tempfile
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def system_setup_dir():
    """Absolute path to the systemSetup/ directory."""
    return Path(__file__).resolve().parent.parent


@pytest.fixture(scope="session")
def input_dir(system_setup_dir):
    """Path to input/ PDB directory."""
    return system_setup_dir / "input"


@pytest.fixture(scope="session")
def pdb_9jjh(input_dir):
    """Path to the smallest example PDB (9jjh)."""
    p = input_dir / "9jjh.pdb"
    if not p.exists():
        pytest.skip("9jjh.pdb not found")
    return str(p)


@pytest.fixture(scope="session")
def pdb_1k3v(input_dir):
    """Path to the primary test PDB (1k3v - PPV)."""
    p = input_dir / "1k3v.pdb"
    if not p.exists():
        pytest.skip("1k3v.pdb not found")
    return str(p)


@pytest.fixture
def workspace():
    """Temporary workspace directory that is cleaned up after each test."""
    with tempfile.TemporaryDirectory(prefix="capsacin_test_") as tmp:
        yield tmp
