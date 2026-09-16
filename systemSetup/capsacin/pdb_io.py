"""Read final sliced PDBs, including older capSACIN serial overflow output."""
from io import StringIO
from pathlib import Path

import MDAnalysis as md


def load_sliced_pdb(path):
    """Repair the old eight-column serial block in memory before parsing.

    Older capSACIN output used five columns for ATOM and eight for its serial
    plus spacing. Above 99,999 the serial overwrote the record's sixth column.
    All later columns stayed in place. Preserve those columns and atom order;
    the source file is never rewritten. Like MDAnalysis's PDB writer, keep the
    final five decimal digits of serials that exceed the PDB field width.
    """
    lines = []
    with Path(path).open() as source:
        for line in source:
            if line.startswith("ATOM "):
                serial = line[5:13].strip()
                if 6 <= len(serial) <= 8 and serial.isdecimal():
                    line = f"ATOM  {int(serial) % 100000:5d}  " + line[13:]
            lines.append(line)
    return md.Universe(StringIO("".join(lines)), format="PDB")
