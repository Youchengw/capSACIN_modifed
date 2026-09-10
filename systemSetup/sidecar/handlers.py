"""Sidecar operation handlers.

Each handler receives a params dict, a progress callback, and a cancel event,
and returns a JSON-serializable dict.
"""
from __future__ import annotations

import json
import os
import tempfile
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import MDAnalysis as md
import numpy as np

from capsacin.pipeline import CapsidPipeline
from capsacin.protocol import SliceRequest, LoadResult, PDBInfoResult, PipelineCancelledError
from capsacin import findSymmetryAxes

from .mmcif_writer import write_aligned_mmcif, write_sliced_mmcif


def handle_inspect_structure(
    params: dict,
    progress: Callable[[str, float, str], None] | None = None,
    cancel: threading.Event | None = None,
) -> dict:
    """Read-only structural metadata.

    Input params:
        input_path: str  — absolute path to PDB file

    Returns:
        chains: list of {id, residue_count, min_resid, max_resid}
        n_models: int
        total_atoms: int
        total_residues: int
        valid_residue_ranges: dict[chain_id, [min, max]]
    """
    input_path = params["input_path"]
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"PDB not found: {input_path}")

    if progress:
        progress("inspect", 0.0, "Opening PDB...")

    u = md.Universe(input_path, input_path)
    return _inspect_universe(u, input_path, progress)


def _inspect_universe(u, input_path, progress=None, totals=None):
    orig_chains = list(np.unique(u.select_atoms("protein").chainIDs))
    n_models = u.trajectory.n_frames

    if progress:
        progress("inspect", 0.3, f"Found {n_models} frames, {len(orig_chains)} chains")

    # Per-chain residue ranges (from frame 0)
    chain_info = []
    valid_ranges = {}
    total_atoms = 0
    total_residues = 0

    for chain in orig_chains:
        sel = u.select_atoms(f"protein and chainid {chain}")
        resids = sel.resids
        if len(resids) > 0:
            min_r = int(np.min(resids))
            max_r = int(np.max(resids))
            n_res = len(np.unique(resids))
        else:
            min_r, max_r, n_res = 0, 0, 0

        chain_info.append({
            "id": chain,
            "residue_count": n_res,
            "min_resid": min_r,
            "max_resid": max_r,
            "atom_count": len(sel),
        })
        valid_ranges[chain] = [min_r, max_r]

    # Count total across all frames
    if totals is None:
        for ts in u.trajectory:
            for chain in orig_chains:
                sel = u.select_atoms(f"protein and chainid {chain}")
                total_atoms += len(sel)
                total_residues += len(np.unique(sel.resids))
    else:
        total_atoms, total_residues = totals

    if progress:
        progress("inspect", 1.0, "Inspection complete")

    return {
        "chains": chain_info,
        "n_models": n_models,
        "n_chains_per_frame": len(orig_chains),
        "total_atoms": total_atoms,
        "total_residues": total_residues,
        "valid_residue_ranges": valid_ranges,
        "input_path": input_path,
    }


@dataclass
class _PreviewSource:
    # Shared only within one request; folds run sequentially because the
    # Universe has a mutable current frame. No structures survive the request.
    loaded: LoadResult
    axis_context: tuple
    pdb_info: PDBInfoResult


def _load_preview_source(input_path, progress=None, cancel=None):
    pipeline = CapsidPipeline(cancel_event=cancel)
    if progress:
        progress("load", 0.0, "Loading PDB...")
    loaded = pipeline.load_structure(input_path)
    if progress:
        progress("detect_axis", 0.07, "Finding 2-, 3- and 5-fold axes...")
    pipeline._check_cancelled()
    context = findSymmetryAxes._compute_axis_context(loaded.universe)
    if progress:
        progress("extract_pdb", 0.15, "Preparing shared atom data...")
    pdb_info = pipeline.extract_pdb_info()
    return _PreviewSource(loaded, context, pdb_info)


def handle_prepare_all_previews(params, workspace_base=None, progress=None, cancel=None):
    """Load once and prepare each fold, retaining individual fold failures."""
    if not params.get("auto", True):
        raise ValueError("Preparing all folds requires automatic axis selection.")
    source = _load_preview_source(params["input_path"], progress, cancel)
    inspection = _inspect_universe(
        source.loaded.universe, params["input_path"],
        totals=(source.loaded.total_atoms, source.loaded.total_residues),
    )
    selected = params.get("symmetry", 5)
    folds = [selected] + [fold for fold in (2, 3, 5) if fold != selected]
    if selected not in (2, 3, 5):
        raise ValueError(f"Unsupported symmetry: {selected}")
    previews, errors = {}, {}
    indices = params.get("axis_indices", {})
    for index, fold in enumerate(folds):
        if cancel is not None and cancel.is_set():
            raise PipelineCancelledError("Pipeline cancelled by user.")
        fraction_done = 0.0

        def fold_progress(stage, fraction, message=""):
            nonlocal fraction_done
            if stage == "prepare":
                fraction_done = max(fraction_done, fraction)
            if progress:
                progress(
                    f"Preparing {fold}-fold ({index + 1}/3)",
                    min(1.0, 0.2 + 0.8 * (index + fraction_done) / 3),
                    message,
                )

        default_index = params.get("axis_index", 0) if fold == selected else 0
        fold_params = dict(params, symmetry=fold, axis_index=indices.get(str(fold), default_index))
        try:
            previews[str(fold)] = handle_prepare_preview(
                fold_params, workspace_base, fold_progress, cancel,
                prepared_source=source,
            )
        except PipelineCancelledError:
            raise
        except Exception as exc:
            errors[str(fold)] = str(exc)
    if progress:
        progress("Previews prepared", 1.0, "All folds processed")
    return {"inspection": inspection, "previews": previews, "errors": errors}


def handle_prepare_preview(
    params: dict,
    workspace_base: Path | None = None,
    progress: Callable[[str, float, str], None] | None = None,
    cancel: threading.Event | None = None,
    *, prepared_source: _PreviewSource | None = None,
) -> dict:
    """Load structure, detect axis, align, and generate viewer mmCIF.

    Input params:
        input_path, symmetry, auto, axis_index, roi_selection,
        roi_frame, ref_indices, legacy_plane_heuristic

    Returns:
        axis_candidates, selected_axis, aligned_cif_path,
        diagnostics, warnings
    """
    input_path = params["input_path"]
    symmetry = params.get("symmetry", 5)
    auto_mode = params.get("auto", True)
    axis_index = params.get("axis_index", 0)
    roi_selection = params.get("roi_selection")
    roi_frame = params.get("roi_frame", 0)
    ref_indices = params.get("ref_indices")
    legacy_plane = params.get("legacy_plane_heuristic", False)

    # Create workspace
    workspace = _get_workspace(workspace_base, "prepare")
    if progress:
        progress("prepare", 0.0, "Creating workspace...")

    source = prepared_source or _load_preview_source(input_path, progress, cancel)
    u = source.loaded.universe
    axis_rows = findSymmetryAxes.list_axes(
        u, symmetry,
        roi_selection=roi_selection,
        roi_frame=roi_frame,
        axis_context=source.axis_context,
        include_references=False,
    )

    candidates = []
    for row in axis_rows:
        candidates.append({
            "axis_index": row["axis_index"],
            "axis": [round(float(v), 6) for v in row["axis"]],
            "score": round(float(row["score"]), 6),
            **{key: row[key] for key in ("ref_index", "ref_frame", "chain_id")
               if row.get(key) is not None},
        })
        if roi_selection:
            candidates[-1]["roi_score"] = round(float(row.get("roi_score", 0)), 6)
            candidates[-1]["roi_angle_deg"] = round(float(row.get("roi_angle_deg", 0)), 3)

    # Run pipeline through alignment
    if progress:
        progress("prepare", 0.1, "Running alignment pipeline...")

    request = SliceRequest(
        input_path=input_path,
        workspace_path=str(workspace),
        symmetry=symmetry,
        weight=0.5,  # dummy — not used for slicing here
        auto=auto_mode,
        axis_index=axis_index,
        roi_selection=roi_selection,
        roi_frame=roi_frame,
        ref_indices=list(ref_indices) if ref_indices else None,
        legacy_plane_heuristic=legacy_plane,
    )

    pipeline = CapsidPipeline(
        progress_callback=progress if progress else None,
        cancel_event=cancel,
        loaded_structure=source.loaded,
    )

    # Run stages 1-5 (load through align)
    load_result = source.loaded
    axis_result = pipeline.detect_axis(
        symmetry=symmetry, auto=auto_mode,
        axis_index=axis_index, roi_selection=roi_selection,
        roi_frame=roi_frame, ref_indices=request.ref_indices,
        legacy_plane_heuristic=legacy_plane,
        axis_context=source.axis_context,
    )
    plane_result = pipeline.select_plane_points(
        axis_result=axis_result, symmetry=symmetry,
        auto_mode=auto_mode, legacy_plane_heuristic=legacy_plane,
    )
    pdb_info = source.pdb_info
    align_result = pipeline.align_coordinates(
        pdb_info=pdb_info, plane_result=plane_result,
        symmetry=symmetry, auto_mode=auto_mode,
        legacy_plane_heuristic=legacy_plane,
    )
    formatted = pipeline.format_aligned_pdb(
        pdb_info=pdb_info, align_result=align_result,
    )

    # Generate viewer mmCIF
    if progress:
        progress("prepare", 0.8, "Generating viewer mmCIF...")

    aligned_cif_path = str(workspace / "aligned_preview.cif")
    write_aligned_mmcif(
        formatted_chains=formatted.chains,
        orig_chains=load_result.orig_chains,
        n_frames=u.trajectory.n_frames,
        output_path=aligned_cif_path,
    )

    # A visual ROI represents the complete local n-fold feature.  Keep the
    # user's roi_frame as the seed and find its symmetry-related MODEL copies
    # around the selected axis for the Mol* viewer.
    roi_symmetry_frames = findSymmetryAxes.find_symmetry_related_frames(
        u,
        plane_result.normal_vector,
        symmetry,
        # A user frame is meaningful only when that ROI actually selected the
        # axis. For global ranking, derive the local copies from the axis.
        seed_frame=roi_frame if roi_selection else None,
        axis_context=source.axis_context,
    )
    pipeline._check_cancelled()

    if progress:
        progress("prepare", 1.0, "Preview ready")

    return {
        "axis_candidates": candidates,
        "selected_axis_index": axis_result.selected_index,
        "selected_axis": [round(float(v), 6) for v in axis_result.axis_dir],
        "aligned_cif_path": aligned_cif_path,
        "workspace_path": str(workspace),
        "diagnostics": {
            "plane_rmsd": plane_result.plane_rmsd,
            "non_collinearity": plane_result.non_collinearity,
            "local_vs_global_angle_deg": plane_result.local_vs_global_angle_deg,
            "bounds": {
                "min": [
                    float(formatted.aligned_coords_df[axis].min())
                    for axis in ("x", "y", "z")
                ],
                "max": [
                    float(formatted.aligned_coords_df[axis].max())
                    for axis in ("x", "y", "z")
                ],
            },
            "roi_symmetry_frames": roi_symmetry_frames,
        },
        "warnings": pipeline._log.copy(),
    }


def handle_run_slice(
    params: dict,
    workspace_base: Path | None = None,
    progress: Callable[[str, float, str], None] | None = None,
    cancel: threading.Event | None = None,
) -> dict:
    """Run the full pipeline and return results.

    Input params:
        input_path, symmetry, weight, auto, axis_index,
        roi_selection, roi_frame, ref_indices, legacy_plane_heuristic

    Returns:
        output_pdb_path, sliced_cif_path, statistics, diagnostics
    """
    input_path = params["input_path"]
    symmetry = params.get("symmetry", 5)
    weight = params.get("weight", 0.5)
    auto_mode = params.get("auto", True)
    axis_index = params.get("axis_index", 0)
    roi_selection = params.get("roi_selection")
    roi_frame = params.get("roi_frame", 0)
    ref_indices = params.get("ref_indices")
    legacy_plane = params.get("legacy_plane_heuristic", False)

    workspace = _get_workspace(workspace_base, "slice")
    if progress:
        progress("slice", 0.0, "Starting full pipeline...")

    request = SliceRequest(
        input_path=input_path,
        workspace_path=str(workspace),
        symmetry=symmetry,
        weight=weight,
        auto=auto_mode,
        axis_index=axis_index,
        roi_selection=roi_selection,
        roi_frame=roi_frame,
        ref_indices=list(ref_indices) if ref_indices else None,
        legacy_plane_heuristic=legacy_plane,
    )

    pipeline = CapsidPipeline(
        progress_callback=progress if progress else None,
        cancel_event=cancel,
    )

    result = pipeline.run_full(request)

    # Generate sliced viewer mmCIF
    sliced_cif_path = str(workspace / "sliced_preview.cif")
    try:
        write_sliced_mmcif(
            cleaned_chains=result.viewer_chains,
            output_path=sliced_cif_path,
            orig_chains=pipeline._orig_chains,
        )
    except Exception as exc:
        pipeline._log_msg(f"[viewer] Could not generate sliced mmCIF: {exc}")
        sliced_cif_path = ""

    return {
        "output_pdb_path": result.output_pdb_path,
        "sliced_cif_path": sliced_cif_path,
        "workspace_path": str(workspace),
        "statistics": {
            "input_atoms": result.input_atoms,
            "output_atoms": result.output_atoms,
            "input_residues": result.input_residues,
            "output_residues": result.output_residues,
            "input_chains": result.input_chains,
            "output_chains": result.output_chains,
            "retained_percentage": round(result.retained_percentage, 2),
        },
        "axis": [round(float(v), 6) for v in result.axis_dir],
        "diagnostics": result.diagnostics,
        "warnings": result.log_lines,
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _get_workspace(base: Path | None, prefix: str) -> Path:
    """Create a unique workspace directory."""
    if base is None:
        base = Path(tempfile.gettempdir()) / "capsacin"
    base.mkdir(parents=True, exist_ok=True)
    # Handler threads and their object IDs can be reused across requests. A
    # real unique directory prevents a later Preview from overwriting a path
    # that the Mol* file cache still associates with an earlier structure.
    return Path(tempfile.mkdtemp(prefix=f"{prefix}_", dir=str(base)))
