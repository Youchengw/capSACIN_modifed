#!/usr/bin/env python3
"""capSACIN: Capsid Surface Abstraction and Computationally-Induced Nanofragmentation.

Steps 1 & 2: Symmetry-based alignment and surface abstraction (slicing).

Usage:
    python sliceCapsid.py --pdb 1k3v --symmetry 5 --weight 0.7 --auto
    python sliceCapsid.py --pdb 1k3v --refindex 959 --symmetry 5 --weight 0.7

The monolithic main() has been refactored to delegate to CapsidPipeline
in capsacin.pipeline.  All existing CLI arguments and output formats are preserved.
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd

import MDAnalysis as md

from capsacin import findSymmetryAxes
from capsacin.pipeline import CapsidPipeline, _print_axis_table
from capsacin.protocol import SliceRequest

# Ensure output directory exists (legacy behaviour)
if not os.path.exists("output"):
    os.makedirs("output")

# Suppress pandas SettingWithCopyWarning safely
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)


# ---------------------------------------------------------------------------
# Legacy helpers retained at module level for any external consumers
# ---------------------------------------------------------------------------


def split_dataframe(df, chunk_size):
    """Split a DataFrame into chunks of *chunk_size* rows."""
    chunks = list()
    num_chunks = len(df) // chunk_size + 1
    for i in range(num_chunks):
        chunks.append(df[i * chunk_size:(i + 1) * chunk_size])
    return chunks


def rotationMatrix(vec1, vec2):
    """Rodrigues rotation matrix from vec1 to vec2.

    Delegates to capsacin.pipeline.rotation_matrix for the actual
    computation; kept here for backward compatibility.
    """
    from capsacin.pipeline import rotation_matrix
    return rotation_matrix(np.asarray(vec1), np.asarray(vec2))


def decode_artificial_chain(artificial_chain, orig_chains):
    """Map a flattened CapSACIN chain ID back to (frame, original chain, flat_index)."""
    from capsacin.pipeline import decode_artificial_chain as _dec
    return _dec(artificial_chain, orig_chains)


# ---------------------------------------------------------------------------
# Main CLI entry point
# ---------------------------------------------------------------------------


def main(args):
    """Run the capSACIN pipeline with CLI arguments.

    Delegates to CapsidPipeline.run_cli() for all computation.
    Legacy stdout formatting is preserved.
    """
    symmetry = args.symmetry
    if symmetry not in (2, 3, 5):
        raise ValueError(
            f"Unsupported symmetry {symmetry}; expected one of 2, 3, or 5."
        )

    output_title = args.pdb
    auto_mode = args.auto
    axis_index = args.axis_index
    roi_selection = args.roi_selection
    roi_frame = args.roi_frame
    weight = args.weight
    legacy_plane = args.legacy_plane_heuristic
    indices_vmd = args.refindex
    plot = args.plot

    pdb_path = f"input/{output_title}.pdb"

    # --- list-axes without --auto: print and exit ---
    if args.list_axes and not auto_mode:
        u = md.Universe(pdb_path, pdb_path)
        axis_rows = findSymmetryAxes.list_axes(
            u, symmetry,
            roi_selection=roi_selection,
            roi_frame=roi_frame,
        )
        _print_axis_table(axis_rows, has_roi=bool(roi_selection))
        return

    if roi_selection and not auto_mode:
        print("[roi] --roi-selection is only used with --auto or --list-axes.")

    # Build request and run through CapsidPipeline
    request = SliceRequest(
        input_path=pdb_path,
        workspace_path="output",
        symmetry=symmetry,
        weight=weight,
        auto=auto_mode,
        axis_index=axis_index,
        roi_selection=roi_selection,
        roi_frame=roi_frame,
        ref_indices=list(indices_vmd) if indices_vmd else None,
        legacy_plane_heuristic=legacy_plane,
        plot=plot,
    )

    pipeline = CapsidPipeline()
    result = pipeline.run_full(request)

    # Print log messages (preserve legacy stdout behaviour)
    for msg in result.log_lines:
        print(msg)

    print(
        f"Retained: {result.retained_percentage:.1f}% "
        f"({result.output_atoms}/{result.input_atoms} atoms, "
        f"{result.output_residues} residues, "
        f"{result.output_chains} chains)"
    )


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="capSACIN: symmetry alignment and surface abstraction"
    )
    parser.add_argument(
        "--pdb", type=str, default="9jjh",
        help="PDB filename prefix (file read from input/{pdb}.pdb)",
    )
    parser.add_argument(
        "--symmetry", "--sym", type=int, choices=(2, 3, 5), default=5,
        help="Symmetry fold (2, 3, or 5)",
    )
    parser.add_argument(
        "--plot", action="store_true",
        help="Enable matplotlib 3D diagnostic plots",
    )
    parser.add_argument(
        "--refindex", type=int, nargs="+", default=[1058],
        help="0-based reference atom indices for manual plane definition",
    )
    parser.add_argument(
        "--weight", type=float, default=0.5,
        help="Slicing weight from 0 (keep all) to 1 (remove all)",
    )
    parser.add_argument(
        "--auto", action="store_true",
        help="Automatically detect symmetry axis and reference atoms",
    )
    parser.add_argument(
        "--axis-index", type=int, default=0,
        help="Which candidate axis to use (0=best; ROI-ranked when "
             "--roi-selection is set)",
    )
    parser.add_argument(
        "--roi-selection", type=str, default=None,
        help="MDAnalysis selection string for ROI-aware axis ranking",
    )
    parser.add_argument(
        "--roi-frame", type=int, default=0,
        help="Trajectory frame (MODEL) containing the ROI copy",
    )
    parser.add_argument(
        "--list-axes", action="store_true",
        help="Print candidate axes table and exit (unless --auto is also set)",
    )
    parser.add_argument(
        "--legacy-plane-heuristic", action="store_true",
        help="Use original 3-point cross-product + z-perturbation "
             "plane-fitting method (backward compatible)",
    )

    args = parser.parse_args()
    main(args)
