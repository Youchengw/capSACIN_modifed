"""Numerical boundaries, PDB/ITP indexing and viewer/export correspondence."""
import threading
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pytest

from capsacin.protocol import PipelineCancelledError
from capsacin.restraints import calculate_restraints
from sidecar.handlers import handle_generate_restraints
from sidecar.main import SidecarServer


def sliced_pdb(tmp_path, legacy_overflow=False):
    path = tmp_path / "slice with spaces.pdb"
    rows = []
    for i, (chain, z) in enumerate([
        ("A", 0), ("A", 30), ("A", 15), ("A", 29.999),
        ("B", 31), ("B", 22.5), ("C", 40),
    ], 1):
        rows.append(f"ATOM  {i:5d}  CA  ALA {chain}{i:4d}    "
                    f"{0:8.3f}{0:8.3f}{z:8.3f}{1:6.2f}{25:6.2f}           C  ")
    if legacy_overflow:
        rows = [f"ATOM {99998 + i:6d}  " + row[13:] for i, row in enumerate(rows)]
    path.write_text("\n".join(rows) + "\nEND\n")
    return path


def test_default_boundaries_and_original_sigmoid():
    z = np.array([0, 14.999, 15, 18.75, 22.5, 26.25, 29.999, 30, 45])
    regions, k = calculate_restraints(z, strength=0.5)
    assert regions.tolist() == [1, 1, 2, 2, 2, 2, 2, 3, 3]
    assert k[:2].tolist() == [500, 500]
    transition = z[2:7]
    normalized = 2 * (transition - transition.min()) / (transition.max() - transition.min()) - 1
    expected = np.round((1 - 1 / (1 + np.exp(-normalized / 0.1))) * 500, 3)
    np.testing.assert_array_equal(k[2:7], expected)
    assert np.all(np.diff(k) <= 0)
    assert k[-2:].tolist() == [0, 0]


def test_empty_and_degenerate_transition():
    regions, k = calculate_restraints([0, 10, 30])
    assert regions.tolist() == [1, 1, 3]
    assert k.tolist() == [1000, 1000, 0]
    with pytest.raises(ValueError, match="one distinct"):
        calculate_restraints([0, 20, 20, 30])


@pytest.mark.parametrize("settings", [
    {"full_until": -1}, {"free_from": 15}, {"full_until": 40},
    {"strength": -1}, {"strength": float("nan")}, {"free_from": float("inf")},
])
def test_invalid_settings(settings):
    with pytest.raises(ValueError, match="finite values"):
        calculate_restraints([0, 18, 25, 35], **settings)


@pytest.mark.parametrize("legacy_overflow", [False, True])
def test_exports_match_viewer_and_per_chain_source_order(tmp_path, monkeypatch, legacy_overflow):
    source = sliced_pdb(tmp_path, legacy_overflow)
    before = source.read_bytes()
    monkeypatch.setenv("CAPSACIN_WORKSPACE", str(tmp_path / "workspaces"))
    server = SidecarServer()
    result = server._dispatch("restraints", "generate_restraints",
                              {"input_path": str(source)}, threading.Event())
    assert source.read_bytes() == before
    assert result["region_counts"] == {"full": 1, "transition": 3, "free": 3}
    assert result["settings"]["strength"] == 1
    assert result["z_range"] == [0, 40]
    files = {f["chain_id"]: f for f in result["files"]}
    a = np.loadtxt(files["A"]["path"], skiprows=2)
    b = np.loadtxt(files["B"]["path"], skiprows=2, ndmin=2)
    assert a[:, 0].tolist() == [1, 3, 4]
    assert b[:, 0].tolist() == [2]
    assert a[0].tolist() == [1, 1, 1000, 1000, 1000]
    assert len(Path(files["C"]["path"]).read_text().splitlines()) == 2
    assert files["C"]["restrained_atoms"] == 0
    viewer_rows = [row.split() for row in Path(result["viewer_cif_path"]).read_text().splitlines()
                   if row.startswith("ATOM ")]
    assert [int(row[-1]) for row in viewer_rows] == [1, 3, 2, 2, 3, 2, 3]
    assert [row[5] for row in viewer_rows] == ["A", "A", "A", "A", "B", "B", "C"]
    with ZipFile(result["archive_path"]) as archive:
        assert len(archive.namelist()) == 4
        for f in files.values():
            assert archive.read(f["filename"]) == Path(f["path"]).read_bytes()
        assert "1-based within each chain" in archive.read("README.txt").decode()


def test_custom_boundaries_and_cancelled_generation(tmp_path):
    source = sliced_pdb(tmp_path)
    first = handle_generate_restraints({"input_path": str(source)}, tmp_path)
    second = handle_generate_restraints({"input_path": str(source), "full_until": 20,
                                          "free_from": 35, "strength": 2}, tmp_path)
    assert second["region_counts"] == {"full": 2, "transition": 4, "free": 1}
    assert first["archive_path"] != second["archive_path"]
    assert Path(first["archive_path"]).is_file()
    event = threading.Event()
    event.set()
    with pytest.raises(PipelineCancelledError):
        handle_generate_restraints({"input_path": str(source)}, tmp_path, cancel=event)
