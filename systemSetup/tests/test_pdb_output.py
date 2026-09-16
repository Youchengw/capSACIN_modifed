"""PDB serial overflow must not hide atoms or collide in the mmCIF viewer."""
import MDAnalysis as md
import numpy as np
import pandas as pd

from capsacin.formatPDB import formatPDB
from sidecar.mmcif_writer import write_aligned_mmcif, write_sliced_mmcif


def test_pdb_serial_boundary_preserves_coordinates_and_viewer_ids(tmp_path):
    serials = [1, 99999, 100000, 100001, 129180]
    coords = np.array([[1.125, 2.25, i * 10.0] for i in range(len(serials))])
    rows = pd.DataFrame([
        {"atom": "ATOM", "idx": serial, "name": "CA", "resname": "ALA",
         "chain": "A" if i < 3 else "B", "resids": i + 1,
         "x": coords[i, 0], "y": coords[i, 1], "z": coords[i, 2],
         "occ": 1.0, "b": 25.0, "type": "C"}
        for i, serial in enumerate(serials)
    ])
    formatted = formatPDB(rows)
    path = tmp_path / "large-serials.pdb"
    np.savetxt(path, formatted, fmt="%s", delimiter="")
    lines = path.read_text().splitlines()
    assert all(line[:6] == "ATOM  " for line in lines)
    assert [int(line[6:11]) for line in lines] == [1, 99999, 0, 1, 29180]
    u = md.Universe(str(path))
    assert len(u.atoms) == 5
    np.testing.assert_array_equal(u.atoms.positions, coords)
    assert u.atoms.chainIDs.tolist() == ["A", "A", "A", "B", "B"]
    assert u.atoms.resids.tolist() == [1, 2, 3, 4, 5]

    # PDB serial 1 appears twice, but mmCIF IDs must remain unique.
    aligned = tmp_path / "aligned.cif"
    sliced = tmp_path / "sliced.cif"
    write_aligned_mmcif(formatted, ["A", "B"], 1, str(aligned))
    write_sliced_mmcif(formatted, str(sliced), ["A", "B"])
    for cif in [aligned, sliced]:
        ids = [line.split()[1] for line in cif.read_text().splitlines()
               if line.startswith("ATOM ")]
        assert ids == ["1", "2", "3", "4", "5"]
