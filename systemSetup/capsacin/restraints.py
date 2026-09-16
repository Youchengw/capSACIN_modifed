"""Position restraints shared by the CLI and Studio.

Heights are rounded z coordinates in the final sliced PDB. The transition uses
the observed minimum/maximum z in that region, as in the original CLI.
"""
from pathlib import Path

import numpy as np

from .pdb_io import load_sliced_pdb


REGION_NAMES = ("full", "transition", "free")


def calculate_restraints(z, full_until=15.0, free_from=30.0, strength=1.0):
    settings = np.asarray([full_until, free_from, strength], dtype=float)
    if not np.isfinite(settings).all() or not 0 <= full_until < free_from or strength < 0:
        raise ValueError("Use finite values: 0 ≤ full boundary < free boundary and strength ≥ 0.")
    z = np.asarray(z, dtype=float)
    if z.ndim != 1 or not len(z) or not np.isfinite(z).all():
        raise ValueError("Restraints require a non-empty slice with finite atom heights.")
    regions = np.where(z < full_until, 1, np.where(z < free_from, 2, 3))
    factors = np.zeros(len(z))
    factors[regions == 1] = 1.0
    transition = z[regions == 2]
    if len(transition):
        span = np.ptp(transition)
        if span == 0:
            raise ValueError("The transition region has only one distinct atom height; "
                             "adjust its boundaries to define the sigmoid normalization.")
        normalized = (transition - transition.min()) * 2.0 / span - 1.0
        factors[regions == 2] = 1.0 - 1.0 / (1.0 + np.exp(-normalized / 0.1))
    return regions, np.round(factors * 1000.0 * strength, 3)


def generate_restraint_files(input_path, output_dir, *, full_until=15.0,
                             free_from=30.0, strength=1.0, prefix=None,
                             check_cancelled=lambda: None):
    check_cancelled()
    u = load_sliced_pdb(input_path)
    if len(u.trajectory) != 1:
        raise ValueError("Generate restraints from the final single-model sliced PDB.")
    atoms = u.atoms
    # Match the CLI's MDAnalysis float32 rounding before normalization.
    z = np.round(atoms.positions[:, 2], 3).astype(float)
    regions, constants = calculate_restraints(z, full_until, free_from, strength)
    chain_ids = atoms.chainIDs
    if any(not str(chain).strip() for chain in chain_ids):
        raise ValueError("The sliced PDB must have a chain ID for every atom.")
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    prefix = prefix or Path(input_path).stem
    files = []
    for chain in np.unique(chain_ids):
        check_cancelled()
        members = np.flatnonzero(chain_ids == chain)
        included = regions[members] != 3
        indices = np.arange(1, len(members) + 1)[included]
        k = constants[members][included]
        # PDB chain IDs can be punctuation; encode unsafe filename characters.
        safe_chain = str(chain) if str(chain).isalnum() else f"U{ord(chain):04X}"
        filename = f"{prefix}-posre-scaled_Protein_chain_{safe_chain}.itp"
        path = directory / filename
        rows = np.column_stack((indices, np.ones(len(indices)), k, k, k))
        np.savetxt(path, rows, fmt=["%d", "%d", "%.3f", "%.3f", "%.3f"],
                   header="[ position_restraints ]\n; atom  type      fx      fy      fz",
                   comments="")
        files.append({"chain_id": str(chain), "path": str(path), "filename": filename,
                      "atom_count": len(members), "restrained_atoms": len(indices)})
    return {"regions": regions, "force_constants": constants, "z": z, "files": files}
