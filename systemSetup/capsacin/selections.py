"""Shared MDAnalysis selections used by the capSACIN pipeline.

MDAnalysis' built-in ``protein`` selection recognizes standard amino-acid
names, but it does not classify CHARMM's protonated glutamate name ``GLUP`` as
protein.  capSACIN must keep that residue in the same atom set at every stage;
otherwise it disappears before the geometric slice is evaluated.
"""

from __future__ import annotations


ADDITIONAL_PROTEIN_RESNAMES = ("GLUP",)

_additional_resnames = " ".join(ADDITIONAL_PROTEIN_RESNAMES)
PROTEIN_SELECTION = f"(protein or resname {_additional_resnames})"


def protein_selection(*clauses: str) -> str:
    """Return the capSACIN protein selection with optional AND clauses."""
    selection = PROTEIN_SELECTION
    for clause in clauses:
        selection = f"({selection} and ({clause}))"
    return selection
