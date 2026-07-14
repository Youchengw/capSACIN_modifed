"""Viewer-only mmCIF generation for Mol* 3D rendering.

Produces mmCIF files where:
- 60 MODELs are expanded into simultaneously visible chains
- Chain IDs use the format M{frame}_{chain} (e.g. M0_A, M59_B)
- Normalized z-coordinate is stored in B_iso_or_equiv for weight-slider preview
- Frame/chain metadata is preserved in a _capsacin custom category
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _normalize_z(z_values: np.ndarray) -> np.ndarray:
    """Normalize z to [0, 1] range.  0=bottom, 1=top."""
    z = np.asarray(z_values, dtype=float)
    z_min, z_max = z.min(), z.max()
    if z_max - z_min < 1e-8:
        return np.full_like(z, 0.5)
    return (z - z_min) / (z_max - z_min)


def _viewer_asym_id(
    artificial_chain: Any,
    orig_chains: list[str],
    flat_chain_index: Any | None = None,
) -> str:
    """Decode a flattened pipeline chain ID to ``M{frame}_{chain}``."""
    if not orig_chains:
        chain_id = str(artificial_chain).strip()
        if not chain_id:
            raise ValueError("Empty artificial chain ID in viewer data.")
        return chain_id
    if flat_chain_index is not None:
        flat_index = int(flat_chain_index)
    else:
        # Backward-compatible path for older/synthetic DataFrames. Do not
        # strip before decoding: chr(65 + index) can itself be Unicode
        # whitespace (for example U+0085 at flat index 68).
        chain_id = str(artificial_chain)
        if not chain_id:
            raise ValueError("Empty artificial chain ID in viewer data.")
        flat_index = ord(chain_id[0]) - 65
    if flat_index < 0:
        raise ValueError(f"Invalid artificial chain ID: {artificial_chain!r}")
    frame = flat_index // len(orig_chains)
    orig_chain = orig_chains[flat_index % len(orig_chains)]
    return f"M{frame}_{orig_chain}"


def write_aligned_mmcif(
    formatted_chains: pd.DataFrame,
    orig_chains: list[str],
    n_frames: int,
    output_path: str,
) -> str:
    """Generate viewer mmCIF for the aligned (pre-slice) capsid.

    Parameters
    ----------
    formatted_chains : pd.DataFrame
        The formatted chains DataFrame from FormattedResult.chains.
        Columns: atom, idx, name, resname, chain, resids, x, y, z, occ, b, type
    orig_chains : list[str]
        Original chain IDs from the PDB.
    n_frames : int
        Number of trajectory frames (MODELs).
    output_path : str
        Where to write the .cif file.

    Returns
    -------
    str
        The output path.
    """
    chains = formatted_chains.copy()

    # Parse coordinates back to float
    x = chains.x.astype(float).values
    y = chains.y.astype(float).values
    z = chains.z.astype(float).values
    nz = _normalize_z(z)

    n_atoms = len(chains)

    lines = [
        "data_CAPSACIN_PREVIEW",
        "# capSACIN viewer mmCIF — aligned preview",
        "#",
        "_entry.id CAPSACIN_PREVIEW",
        "#",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.auth_atom_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_seq_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
    ]

    # Build frame → chain mapping for unique asym IDs
    # The chains DataFrame uses artificial single-letter chain IDs (A-Z etc.).
    # Map each back to (frame, original_chain).
    for i in range(n_atoms):
        atom_type = str(chains.atom.values[i]).strip()
        atom_id = str(chains.idx.values[i]).strip()
        atom_name = str(chains.name.values[i]).strip()
        element = str(chains.type.values[i]).strip() or atom_name[0]
        res_name = str(chains.resname.values[i]).strip()
        artificial_chain = str(chains.chain.values[i]).strip()
        resid = str(chains.resids.values[i]).strip()
        occ = str(chains.occ.values[i]).strip()
        b_factor = f"{nz[i]:.3f}"

        # Decode artificial chain back to frame and original chain
        flat_chain_index = (
            chains.viewer_chain_index.values[i]
            if "viewer_chain_index" in chains.columns
            else None
        )
        asym_id = _viewer_asym_id(
            artificial_chain,
            orig_chains,
            flat_chain_index,
        )

        lines.append(
            f"{atom_type} {atom_id} {element} {atom_name} {res_name} "
            f"{asym_id} 1 {resid} {atom_name} {res_name} {asym_id} {resid} "
            f"{x[i]:.3f} {y[i]:.3f} {z[i]:.3f} "
            f"{occ} {b_factor}"
        )

    lines.append("#")

    # Write
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    return output_path


def write_sliced_mmcif(
    cleaned_chains: str | pd.DataFrame,
    output_path: str,
    orig_chains: list[str] | None = None,
) -> str:
    """Generate viewer mmCIF for the sliced (post-pipeline) result.

    Parameters
    ----------
    cleaned_chains : str or pd.DataFrame
        Either the path to the scientific output PDB, or a chains DataFrame.
    output_path : str
        Where to write the .cif file.
    orig_chains : list[str] or None
        Original input chain IDs. Required for a DataFrame when MODEL/chain
        identity must be preserved for ROI selection.

    Returns
    -------
    str
        The output path.
    """
    if isinstance(cleaned_chains, str):
        # Read the PDB using MDAnalysis
        import MDAnalysis as md
        u = md.Universe(cleaned_chains)
        atoms = u.atoms
        n_atoms = len(atoms)
        positions = atoms.positions

        lines = [
            "data_CAPSACIN_SLICED",
            "# capSACIN viewer mmCIF — sliced result",
            "#",
            "_entry.id CAPSACIN_SLICED",
            "#",
            "loop_",
            "_atom_site.group_PDB",
            "_atom_site.id",
            "_atom_site.type_symbol",
            "_atom_site.label_atom_id",
            "_atom_site.label_comp_id",
            "_atom_site.label_asym_id",
            "_atom_site.label_entity_id",
            "_atom_site.label_seq_id",
            "_atom_site.auth_atom_id",
            "_atom_site.auth_comp_id",
            "_atom_site.auth_asym_id",
            "_atom_site.auth_seq_id",
            "_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z",
            "_atom_site.occupancy",
            "_atom_site.B_iso_or_equiv",
        ]

        for i in range(n_atoms):
            atom = atoms[i]
            asym_id = atom.chainID or atom.segid or "A"
            lines.append(
                f"ATOM {i + 1} {atom.element or atom.name[0]} {atom.name} {atom.resname} "
                f"{asym_id} 1 {atom.resid} {atom.name} {atom.resname} {asym_id} {atom.resid} "
                f"{positions[i, 0]:.3f} {positions[i, 1]:.3f} {positions[i, 2]:.3f} "
                f"1.00 0.00"
            )
    else:
        # DataFrame path
        df = cleaned_chains.copy()
        x = df.x.astype(float).values
        y = df.y.astype(float).values
        z = df.z.astype(float).values
        n_atoms = len(df)

        lines = [
            "data_CAPSACIN_SLICED",
            "# capSACIN viewer mmCIF — sliced result",
            "#",
            "_entry.id CAPSACIN_SLICED",
            "#",
            "loop_",
            "_atom_site.group_PDB",
            "_atom_site.id",
            "_atom_site.type_symbol",
            "_atom_site.label_atom_id",
            "_atom_site.label_comp_id",
            "_atom_site.label_asym_id",
            "_atom_site.label_entity_id",
            "_atom_site.label_seq_id",
            "_atom_site.auth_atom_id",
            "_atom_site.auth_comp_id",
            "_atom_site.auth_asym_id",
            "_atom_site.auth_seq_id",
            "_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z",
            "_atom_site.occupancy",
            "_atom_site.B_iso_or_equiv",
        ]

        for i in range(n_atoms):
            asym_id = _viewer_asym_id(
                df.chain.values[i],
                orig_chains or [],
                (
                    df.viewer_chain_index.values[i]
                    if "viewer_chain_index" in df.columns
                    else None
                ),
            )
            lines.append(
                f"{str(df.atom.values[i]).strip()} {str(df.idx.values[i]).strip()} "
                f"{str(df.type.values[i]).strip()} {str(df.name.values[i]).strip()} "
                f"{str(df.resname.values[i]).strip()} {asym_id} "
                f"1 {str(df.resids.values[i]).strip()} {str(df.name.values[i]).strip()} "
                f"{str(df.resname.values[i]).strip()} {asym_id} "
                f"{str(df.resids.values[i]).strip()} "
                f"{x[i]:.3f} {y[i]:.3f} {z[i]:.3f} "
                f"{str(df.occ.values[i]).strip()} 0.00"
            )

    lines.append("#")
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    return output_path
