#!/usr/bin/env python3
"""Structured, callable capSACIN pipeline.

Extracts the monolithic logic from sliceCapsid.py into a CapsidPipeline class
with discrete stages, progress callbacks, cancellation, and structured results.
"""

from __future__ import annotations

import copy
import os
import threading
import warnings
from typing import Any, Callable

import numpy as np
from numpy import savetxt
import pandas as pd
from tqdm import tqdm

import MDAnalysis as md
from MDAnalysis.analysis import contacts

from capsacin import definePlane, formatPDB, createDictionary, findSymmetryAxes
from capsacin.definePlane import (
    fit_plane_normal,
    select_5fold_ring,
    compute_pentagon_diagnostics,
)
from capsacin.protocol import (
    AlignResult,
    AxisCandidate,
    AxisResult,
    CleanupResult,
    FormattedResult,
    LoadResult,
    PDBInfoResult,
    PipelineCancelledError,
    PlaneResult,
    SliceDataResult,
    SliceOutput,
    SliceRequest,
)

# Suppress pandas SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=pd.errors.SettingWithCopyWarning)


# ---------------------------------------------------------------------------
# Pure helpers (module-level, unchanged from sliceCapsid.py)
# ---------------------------------------------------------------------------


def rotation_matrix(vec1: np.ndarray, vec2: np.ndarray) -> np.ndarray:
    """Rodrigues rotation matrix from vec1 to vec2."""
    a = (vec1 / np.linalg.norm(vec1)).reshape(3)
    b = (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.clip(np.dot(a, b), -1.0, 1.0)
    s = np.linalg.norm(v)
    if s < 1e-12:
        if c > 0:
            return np.eye(3)
        helper = np.array([1.0, 0.0, 0.0])
        if abs(a[0]) > 0.9:
            helper = np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, helper)
        axis = axis / np.linalg.norm(axis)
        kmat = np.array([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0],
        ])
        return np.eye(3) + 2.0 * kmat.dot(kmat)
    kmat = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0],
    ])
    return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))


def decode_artificial_chain(artificial_chain: Any, orig_chains: list) -> tuple:
    """Map a flattened CapSACIN chain ID back to (frame, original_chain, flat_index)."""
    chain_id = str(artificial_chain)
    if len(chain_id) == 0:
        raise ValueError("Empty artificial chain ID.")
    flat_index = ord(chain_id[0]) - 65
    if flat_index < 0:
        raise ValueError(f"Invalid artificial chain ID: {artificial_chain!r}")
    n_chains = len(orig_chains)
    frame = flat_index // n_chains
    chain_slot = flat_index % n_chains
    return frame, orig_chains[chain_slot], flat_index


# ---------------------------------------------------------------------------
# CapsidPipeline
# ---------------------------------------------------------------------------


class CapsidPipeline:
    """Callable pipeline wrapping sliceCapsid.py logic.

    Parameters
    ----------
    progress_callback : callable or None
        Called as ``callback(stage_name: str, fraction: float)`` at stage boundaries.
    cancel_event : threading.Event or None
        When set, the next cancellation checkpoint raises PipelineCancelledError.
    """

    def __init__(
        self,
        progress_callback: Callable[[str, float], None] | None = None,
        cancel_event: threading.Event | None = None,
        loaded_structure: LoadResult | None = None,
    ):
        self._progress = progress_callback
        self._cancel = cancel_event or threading.Event()

        # Internal state populated by stages
        self._u: md.Universe | None = None
        self._orig_chains: list[str] = []
        self._n_monomers: int = 0
        self._load_result: LoadResult | None = None
        self._log: list[str] = []
        if loaded_structure is not None:
            self._u = loaded_structure.universe
            self._orig_chains = loaded_structure.orig_chains
            self._n_monomers = loaded_structure.n_monomers
            self._load_result = loaded_structure

    # ------------------------------------------------------------------
    # Cancellation / progress helpers
    # ------------------------------------------------------------------

    def _check_cancelled(self) -> None:
        if self._cancel.is_set():
            raise PipelineCancelledError("Pipeline cancelled by user.")

    def _emit_progress(self, stage: str, fraction: float) -> None:
        if self._progress:
            self._progress(stage, fraction)

    def _log_msg(self, msg: str) -> None:
        self._log.append(msg)

    # ------------------------------------------------------------------
    # Stage 1 — Load PDB
    # ------------------------------------------------------------------

    def load_structure(self, pdb_path: str) -> LoadResult:
        """Open a PDB file and extract chain metadata.

        Parameters
        ----------
        pdb_path : str
            Absolute or relative path to the multi-MODEL PDB file.

        Returns
        -------
        LoadResult
        """
        self._check_cancelled()
        self._emit_progress("load", 0.0)

        u = md.Universe(pdb_path, pdb_path)
        orig_chains = list(np.unique(u.select_atoms("protein").chainIDs))
        n_monomers = len(orig_chains)

        # Count total atoms and residues
        total_atoms = 0
        total_residues = 0
        for ts in u.trajectory:
            for chain in orig_chains:
                sel = u.select_atoms(f"protein and chainid {chain}")
                total_atoms += len(sel)
                total_residues += len(np.unique(sel.resids))

        self._u = u
        self._orig_chains = orig_chains
        self._n_monomers = n_monomers

        load_result = LoadResult(
            universe=u,
            orig_chains=orig_chains,
            n_monomers=n_monomers,
            total_atoms=total_atoms,
            total_residues=total_residues,
        )
        self._load_result = load_result
        self._emit_progress("load", 1.0)
        return load_result

    # ------------------------------------------------------------------
    # Stage 2 — Detect symmetry axis
    # ------------------------------------------------------------------

    def detect_axis(
        self,
        symmetry: int,
        *,
        auto: bool = True,
        axis_index: int = 0,
        roi_selection: str | None = None,
        roi_frame: int = 0,
        ref_indices: list[int] | None = None,
        legacy_plane_heuristic: bool = False,
        list_axes: bool = False,
        axis_context=None,
    ) -> AxisResult:
        """Detect or validate the symmetry axis.

        In auto mode, calls ``findSymmetryAxes.auto_detect_reference``.
        In manual mode, validates the user-supplied ``ref_indices``.

        Parameters
        ----------
        symmetry : int
            2, 3, or 5.
        auto : bool
            Enable automatic axis + reference detection.
        axis_index : int
            Which candidate axis to use (0 = best; ROI-ranked when roi_selection set).
        roi_selection : str or None
            MDAnalysis atom selection string for ROI-aware ranking.
        roi_frame : int
            Trajectory frame containing the ROI.
        ref_indices : list[int] or None
            Manual reference atom indices (used when auto=False).
        legacy_plane_heuristic : bool
            Legacy plane behaviour flag (passed through to diagnostics).
        list_axes : bool
            If True and auto is False, returns candidate list without running pipeline.

        Returns
        -------
        AxisResult
        """
        self._check_cancelled()
        self._emit_progress("detect_axis", 0.0)

        if self._u is None:
            raise RuntimeError("load_structure() must be called before detect_axis()")

        if symmetry not in (2, 3, 5):
            raise ValueError(f"Unsupported symmetry {symmetry}; expected 2, 3, or 5.")

        candidates: list[AxisCandidate] = []
        axis_dir: list[float] = [0.0, 0.0, 1.0]
        result_ref_indices: list[int] = []
        ref_frames: list[int] = []
        selected_index: int = axis_index

        # --- list_axes mode (no pipeline execution) ---
        if list_axes:
            axis_rows = findSymmetryAxes.list_axes(
                self._u, symmetry,
                roi_selection=roi_selection,
                roi_frame=roi_frame,
            )
            for row in axis_rows:
                candidates.append(AxisCandidate(
                    axis=list(row["axis"]),
                    score=row["score"],
                    axis_index=row["axis_index"],
                    chain_id=row.get("chain_id"),
                    ref_frame=row.get("ref_frame"),
                    ref_index=row.get("ref_index"),
                    roi_score=row.get("roi_score"),
                    roi_angle_deg=row.get("roi_angle_deg"),
                    roi_line_distance=row.get("roi_line_distance"),
                ))
            self._emit_progress("detect_axis", 1.0)
            return AxisResult(
                axis_dir=axis_dir,
                ref_indices=result_ref_indices,
                ref_frames=ref_frames,
                candidates=candidates,
                selected_index=selected_index,
            )

        # --- Auto mode ---
        if auto:
            auto_ref, axis_dir_arr, auto_frames = findSymmetryAxes.auto_detect_reference(
                self._u, symmetry,
                axis_index=axis_index,
                roi_selection=roi_selection,
                roi_frame=roi_frame,
                axis_context=axis_context,
            )
            result_ref_indices = list(auto_ref)
            axis_dir = list(float(v) for v in axis_dir_arr)
            ref_frames = list(auto_frames)
            selected_index = axis_index

            if roi_selection:
                self._log_msg(f"[auto] ROI selection: {roi_selection!r} on frame {roi_frame}")
                self._log_msg(f"[auto] Using ROI-ranked {symmetry}-fold axis at rank {axis_index}")
            self._log_msg(f"[auto] Detected {symmetry}-fold axis: {axis_dir}")
            self._log_msg(f"[auto] Reference atom indices: {result_ref_indices}")
            self._log_msg(f"[auto] Reference frames: {ref_frames}")
        else:
            # Manual mode — use user-supplied ref_indices
            if ref_indices is None or len(ref_indices) == 0:
                raise ValueError("Manual mode requires --refindex to be set.")
            # Only pointA index is used; pointB/C found via distance ranking
            result_ref_indices = list(ref_indices)
            # axis_dir stays [0,0,1]; will be recomputed from plane fit
            ref_frames = [0, 0, 0]

        self._emit_progress("detect_axis", 1.0)
        return AxisResult(
            axis_dir=axis_dir,
            ref_indices=result_ref_indices,
            ref_frames=ref_frames,
            candidates=candidates,
            selected_index=selected_index,
        )

    # ------------------------------------------------------------------
    # Stage 3 — Select plane points
    # ------------------------------------------------------------------

    def select_plane_points(
        self,
        axis_result: AxisResult,
        symmetry: int,
        *,
        auto_mode: bool = True,
        legacy_plane_heuristic: bool = False,
        plot: bool = False,
    ) -> PlaneResult:
        """Select reference points for plane fitting.

        This is the most complex stage due to fold-specific branching.
        In auto mode, uses the detected monomer positions directly.
        In manual mode, uses nearest-neighbour distance ranking.

        Returns
        -------
        PlaneResult
            Contains pointA/B/C, plane_points, normal_vector, and diagnostics.
        """
        self._check_cancelled()
        self._emit_progress("plane_points", 0.0)

        u = self._u
        if u is None:
            raise RuntimeError("load_structure() must be called first.")

        orig_chains = self._orig_chains
        ref_indices = axis_result.ref_indices
        ref_frames = axis_result.ref_frames
        axis_dir = np.array(axis_result.axis_dir, dtype=float)

        # Default values
        pointA: np.ndarray = np.zeros(3)
        pointB: np.ndarray = np.zeros(3)
        pointC: np.ndarray = np.zeros(3)
        plane_points: np.ndarray = np.zeros((3, 3))
        normal_vector: np.ndarray = np.array([0.0, 0.0, 1.0])
        auto_plane_points: np.ndarray | None = None

        plane_rmsd: float | None = None
        non_collinearity: float | None = None
        local_vs_global_angle_deg: float | None = None

        # --- Get pointA from reference atom ---
        u.trajectory[0]  # reset trajectory
        refPos = u.select_atoms(f"protein and index {ref_indices[0]}")
        pointA = refPos.positions[0].copy()

        # ==================================================================
        # AUTO MODE: pointB/C from detected monomers
        # ==================================================================
        if auto_mode and len(ref_indices) >= 3 and len(ref_frames) >= 3:
            # pointB from second detected monomer
            u.trajectory[ref_frames[1]]
            refPosB_md = u.select_atoms(f"protein and index {ref_indices[1]}")
            pointB = refPosB_md.positions[0].copy()

            # pointC from third detected monomer
            u.trajectory[ref_frames[2]]
            refPosC_md = u.select_atoms(f"protein and index {ref_indices[2]}")
            pointC = refPosC_md.positions[0].copy()

            # Return to first monomer for pointA
            u.trajectory[ref_frames[0]]
            refPosA_md = u.select_atoms(f"protein and index {ref_indices[0]}")
            pointA = refPosA_md.positions[0].copy()

            if legacy_plane_heuristic:
                pointB[2] += 0.1
                pointC[2] -= 0.1
            else:
                # Non-legacy auto: compute local plane diagnostics
                if symmetry == 2:
                    pair_points, pair_diag = findSymmetryAxes.select_2fold_pairs(
                        u,
                        int(ref_indices[0]),
                        ref_frame=ref_frames[0],
                        partner_frame=ref_frames[2],
                        ref_chain_id=refPosA_md.chainIDs[0],
                        partner_chain_id=refPosA_md.chainIDs[0],
                    )
                    if pair_points is not None:
                        auto_plane_points = pair_points
                        pointB = pair_points[1]
                        pointC = pair_points[2]
                        self._log_msg(
                            f"[plane] auto fold=2 2nd-atom: "
                            f"index={pair_diag['b_atom_index']} "
                            f"name={pair_diag['b_name']} "
                            f"resid={pair_diag['b_resid']} "
                            f"partner_frame={pair_diag['partner_frame']} "
                            f"partner_chain={pair_diag['partner_chain_id']} "
                            f"pair_dist={pair_diag['pair_distance']:.1f}A "
                            f"candidates={pair_diag['candidate_count']}"
                        )
                    else:
                        self._log_msg(
                            f"[plane] auto fold=2 2nd-atom selection failed: "
                            f"{pair_diag.get('error', 'unknown')}; using 3-point diagnostic"
                        )

                auto_pts = (
                    auto_plane_points
                    if auto_plane_points is not None
                    else np.array([pointA, pointB, pointC])
                )
                auto_normal, auto_rmsd, auto_sv = fit_plane_normal(
                    auto_pts, reference_direction=axis_dir,
                )
                alignment_dot = np.clip(np.abs(np.dot(auto_normal, axis_dir)), -1.0, 1.0)
                alignment_deg = float(np.rad2deg(np.arccos(alignment_dot)))
                ncl = float(auto_sv[1] / auto_sv[0] if auto_sv[0] > 1e-12 else 0.0)
                self._log_msg(
                    f"[plane] auto fold={symmetry} local_rmsd={auto_rmsd:.4f}A "
                    f"non_coll={ncl:.2f} local_vs_global_angle={alignment_deg:.2f}deg"
                )
                if auto_rmsd > 2.0:
                    self._log_msg(
                        f"[plane] WARNING: auto local plane RMSD {auto_rmsd:.2f}A "
                        f"exceeds 2.0A threshold"
                    )
                if ncl < 0.05:
                    self._log_msg(
                        f"[plane] WARNING: near-collinear auto points "
                        f"(non_coll={ncl:.3f}); normal may be unreliable"
                    )
                if alignment_deg > 5.0:
                    self._log_msg(
                        f"[plane] WARNING: local plane vs global axis angle "
                        f"{alignment_deg:.2f}deg exceeds 5deg threshold"
                    )
                plane_rmsd = float(auto_rmsd)
                non_collinearity = ncl
                local_vs_global_angle_deg = alignment_deg

            # In auto mode, the global axis is used for alignment
            normal_vector = axis_dir.copy()

        # ==================================================================
        # MANUAL MODE: nearest-neighbour distance ranking
        # ==================================================================
        else:
            # Need to build chains DataFrame first for candidate matching
            # We'll build a temporary chains-like structure for candidate selection
            u.trajectory[0]  # reset
            refPos = u.select_atoms(f"protein and index {ref_indices[0]}")

            # Build candidate coordinates by extracting matching atoms across all frames
            candidate_atom_list = []
            candidate_idx_list = []
            candidate_chain_list = []
            for ts in u.trajectory:
                for chain in orig_chains:
                    virus = u.select_atoms(f"protein and chainid {chain}")
                    for atom_pos, atom_idx, atom_name, atom_resname, atom_resid in zip(
                        virus.positions, np.arange(len(virus)), virus.names,
                        virus.resnames, virus.resids,
                    ):
                        if (atom_resname == refPos.resnames[0]
                                and atom_resid == refPos.resids[0]
                                and atom_name == refPos.names[0]):
                            # Construct artificial chain ID matching legacy behavior
                            # We need the flattened chain ID used in the chains DataFrame
                            candidate_atom_list.append(atom_pos)
                            candidate_idx_list.append(atom_idx)
                            # Use a placeholder chain — will be matched via distance
                            candidate_chain_list.append(f"{chain}_{ts.frame}")

            if len(candidate_atom_list) < 2:
                raise RuntimeError(
                    f"Fewer than 2 matching atoms found for manual reference; "
                    f"check refindex {ref_indices[0]}"
                )

            candidate_coords = np.array(candidate_atom_list)
            candidate_indices = np.array(candidate_idx_list)

            dist_matrix = contacts.distance_array(pointA, candidate_coords)
            dist_matrix = np.vstack((candidate_indices, dist_matrix)).T
            dist_sorted = np.argsort(dist_matrix[:, 1])

            if legacy_plane_heuristic or symmetry == 3:
                if symmetry == 3:
                    k1, k2 = 1, 2
                elif symmetry == 5:
                    k1, k2 = 1, 3
                elif symmetry == 2:
                    k1, k2 = 0, 1
                else:
                    k1, k2 = 1, 2

                pointB = candidate_coords[dist_sorted[k1]].copy()
                pointB[2] += 0.1
                pointC = candidate_coords[dist_sorted[k2]].copy()
                pointC[2] -= 0.1
                plane_points = np.array([pointA, pointB, pointC])

            elif symmetry == 5:
                ring = select_5fold_ring(candidate_coords, pointA)
                pointB = ring[1].copy()
                pointC = ring[3].copy()
                plane_points = ring[:5].copy()

            elif symmetry == 2:
                ref_idx = int(refPos.indices[0])
                # Determine partner frame from distance ranking
                partner_idx = int(dist_sorted[1])
                partner_coord_idx = int(candidate_indices[partner_idx])
                # Partner's position in the candidate list tells us which frame
                partner_frame = partner_idx // len(orig_chains)
                if 0 <= partner_frame < u.trajectory.n_frames:
                    pair_points, pair_diag = findSymmetryAxes.select_2fold_pairs(
                        u, ref_idx, ref_frame=0, partner_frame=partner_frame,
                        ref_chain_id=refPos.chainIDs[0],
                        partner_chain_id=orig_chains[partner_idx % len(orig_chains)],
                    )
                    if pair_points is not None:
                        plane_points = pair_points.copy()
                        pointB = pair_points[1].copy()
                        pointC = pair_points[2].copy()
                        self._log_msg(
                            f"[plane] fold=2 2nd-atom: "
                            f"index={pair_diag['b_atom_index']} "
                            f"name={pair_diag['b_name']} "
                            f"resid={pair_diag['b_resid']} "
                            f"partner_frame={pair_diag['partner_frame']} "
                            f"partner_chain={pair_diag['partner_chain_id']} "
                            f"pair_dist={pair_diag['pair_distance']:.1f}A "
                            f"candidates={pair_diag['candidate_count']}"
                        )
                    else:
                        self._log_msg(
                            f"[plane] fold=2 2nd-atom selection failed: "
                            f"{pair_diag.get('error', 'unknown')}; falling back to z-perturbation"
                        )
                        k1, k2 = 0, 1
                        pointB = candidate_coords[dist_sorted[k1]].copy()
                        pointB[2] += 0.1
                        pointC = candidate_coords[dist_sorted[k2]].copy()
                        pointC[2] -= 0.1
                        plane_points = np.array([pointA, pointB, pointC])
                else:
                    self._log_msg(
                        "[plane] fold=2 could not determine partner frame; "
                        "falling back to z-perturbation"
                    )
                    k1, k2 = 0, 1
                    pointB = candidate_coords[dist_sorted[k1]].copy()
                    pointB[2] += 0.1
                    pointC = candidate_coords[dist_sorted[k2]].copy()
                    pointC[2] -= 0.1
                    plane_points = np.array([pointA, pointB, pointC])

            # In manual + non-legacy + non-auto mode, normal_vector will be
            # computed in the alignment stage via SVD.

        # Override plane_points if auto_plane_points was set
        if auto_plane_points is not None:
            plane_points = auto_plane_points

        self._emit_progress("plane_points", 1.0)
        return PlaneResult(
            pointA=pointA,
            pointB=pointB,
            pointC=pointC,
            plane_points=plane_points,
            normal_vector=normal_vector,
            plane_rmsd=plane_rmsd,
            non_collinearity=non_collinearity,
            local_vs_global_angle_deg=local_vs_global_angle_deg,
        )

    # ------------------------------------------------------------------
    # Stage 4 — Extract PDB info (chains DataFrame + coordinates)
    # ------------------------------------------------------------------

    def extract_pdb_info(self) -> PDBInfoResult:
        """Build the chains DataFrame and coordinate array from the Universe.

        Iterates over all trajectory frames and chains, collecting atom-level
        metadata and assembling the flat chains DataFrame used by all subsequent
        stages.

        Returns
        -------
        PDBInfoResult
        """
        self._check_cancelled()
        self._emit_progress("extract_pdb", 0.0)

        u = self._u
        orig_chains = self._orig_chains
        if u is None:
            raise RuntimeError("load_structure() must be called first.")

        # Build allCoords (list of per-chain position arrays)
        all_coords = []
        for ts in tqdm(u.trajectory, desc="Extracting coordinates"):
            self._check_cancelled()
            for chain in orig_chains:
                virus_coords = u.select_atoms(f"protein and chainid {chain}").positions
                all_coords.append(virus_coords)

        built_virus = np.vstack(all_coords)
        capsid_com = np.array([
            np.mean(built_virus[:, 0]),
            np.mean(built_virus[:, 1]),
            np.mean(built_virus[:, 2]),
        ])

        # Build amino acid checks dictionary
        amino_acid_checks = createDictionary.createDictionary(u)

        # Extract per-chain atom metadata
        resids_list, resnames_list, bfactors_list = [], [], []
        names_list, types_list, occupancies_list = [], [], []
        indices_list, atoms_list = [], []
        len_chain = []

        for j, chain in enumerate(orig_chains):
            virus = u.select_atoms(f"protein and chainid {chain}")
            resids_list.append(virus.resids)
            resnames_list.append(virus.resnames)
            bfactors_list.append(virus.tempfactors)
            names_list.append(virus.names)
            types_list.append(virus.types)
            occupancies_list.append(virus.occupancies)
            indices_list.append(np.arange(1, len(all_coords[j]) + 1, 1))
            len_chain.append(len(np.unique(resids_list[j])))
            atom_col = []
            for _ in range(len(all_coords[j])):
                atom_col.append("ATOM")
            atoms_list.append(atom_col)

        n_frames = u.trajectory.n_frames
        N = self._n_monomers * n_frames
        chain_id = [chr(65 + i) for i in range(N)]

        # Stack into chains DataFrame
        chains_list = []
        k = 0
        for i in range(len(all_coords)):
            x = all_coords[i][:, 0]
            y = all_coords[i][:, 1]
            z = all_coords[i][:, 2]

            chain_list_local = []
            for _ in range(len(all_coords[k])):
                chain_list_local.append(chain_id[i])

            index = indices_list[k] + (i * (len(all_coords[k])))
            curr_chain = np.vstack((
                atoms_list[k], index, names_list[k], resnames_list[k],
                chain_list_local, resids_list[k], x, y, z,
                occupancies_list[k], bfactors_list[k], types_list[k],
            )).T
            chains_list.append(curr_chain)
            k += 1
            if k == len(orig_chains):
                k = 0

        chains_arr = np.vstack(chains_list)
        chains = pd.DataFrame(chains_arr, columns=[
            "atom", "idx", "name", "resname", "chain",
            "resids", "x", "y", "z", "occ", "b", "type",
        ])
        # Stable identity for viewer mmCIF generation. The legacy scientific
        # pipeline encodes flattened chains as chr(65 + index); with enough
        # MODEL/chain copies this reaches Unicode whitespace characters. Keep
        # an explicit numeric index so the viewer never has to reverse that
        # fragile one-character encoding.
        chains["viewer_chain_index"] = np.concatenate([
            np.full(len(coords), index, dtype=int)
            for index, coords in enumerate(all_coords)
        ])

        coords = np.vstack((
            chains.x.astype(float), chains.y.astype(float), chains.z.astype(float),
        )).T.copy()

        self._emit_progress("extract_pdb", 1.0)
        return PDBInfoResult(
            chains=chains,
            coords=coords,
            capsid_com=capsid_com,
            amino_acid_checks=amino_acid_checks,
        )

    # ------------------------------------------------------------------
    # Stage 5 — Align coordinates
    # ------------------------------------------------------------------

    def align_coordinates(
        self,
        pdb_info: PDBInfoResult,
        plane_result: PlaneResult,
        symmetry: int,
        *,
        auto_mode: bool = True,
        legacy_plane_heuristic: bool = False,
        rotation_check: bool = True,
        plot: bool = False,
    ) -> AlignResult:
        """Compute alignment rotation and apply to coordinates.

        Uses the plane normal (or detected axis in auto mode) to rotate the
        capsid so the symmetry axis aligns with the z-axis.

        Returns
        -------
        AlignResult
        """
        self._check_cancelled()
        self._emit_progress("align", 0.0)

        u = self._u
        if u is None:
            raise RuntimeError("load_structure() must be called first.")

        coords = pdb_info.coords.copy()
        capsid_com = pdb_info.capsid_com
        pointA = plane_result.pointA
        pointB = plane_result.pointB
        pointC = plane_result.pointC
        plane_points = plane_result.plane_points.copy()
        normal_vector = plane_result.normal_vector.copy()
        axis_dir = np.array(plane_result.normal_vector, dtype=float)

        # Determine which path to use for normal computation
        if auto_mode and len(getattr(self, '_auto_ref_indices', []) or []) > 0:
            # Actually check if we're in a true auto mode with valid ref_indices
            pass  # Will be determined below

        # Build selectPoints and center
        if legacy_plane_heuristic or (
            auto_mode and plane_result.normal_vector is not None
            and np.linalg.norm(plane_result.normal_vector) > 0.5
        ):
            # Legacy path: 3-point DataFrame
            points_df = pd.DataFrame(
                np.vstack((pointA, pointB, pointC)),
                columns=['x', 'y', 'z'],
            )
            x, y, z_pts = points_df.x.values, points_df.y.values, points_df.z.values
            select_points = np.vstack((x, y, z_pts)).T
            com_pt = np.array([np.mean(x), np.mean(y), np.mean(z_pts)])
        else:
            # SVD path: use plane_points directly
            select_points = plane_points.copy()
            com_pt = plane_points.mean(axis=0)

        coords_centered = coords - com_pt
        select_centered = select_points - com_pt

        # Compute normal vector for alignment
        use_auto_axis = bool(
            auto_mode
            and np.linalg.norm(axis_dir) > 0.5
            and not legacy_plane_heuristic
        )

        if use_auto_axis:
            normal = axis_dir
            self._log_msg("[auto] Using detected axis as alignment normal (bypassing 3-point method)")
        elif legacy_plane_heuristic:
            x_s, y_s, z_s = select_centered[:, 0], select_centered[:, 1], select_centered[:, 2]
            normal = definePlane.definePlane(x_s, y_s, z_s)
        else:
            # SVD-based plane fit
            ref_dir = pointA - capsid_com
            # Optional global axis cross-check
            global_axis = None
            if symmetry in (2, 5):
                try:
                    _cc = findSymmetryAxes.compute_capsid_center(u)
                    _mcoms = findSymmetryAxes.compute_monomer_coms(u)
                    _axes = findSymmetryAxes.find_symmetry_axes(
                        _mcoms, _cc, fold=symmetry, n_grid=2000,
                    )
                    if len(_axes) > 0:
                        dots_ref = _axes @ (ref_dir / np.linalg.norm(ref_dir))
                        best_i = int(np.argmax(np.abs(dots_ref)))
                        global_axis = _axes[best_i]
                        if np.dot(global_axis, ref_dir) < 0:
                            global_axis = -global_axis
                except Exception as e:
                    self._log_msg(f"[plane] Could not compute global axis: {e}")

            orientation_ref = global_axis if global_axis is not None else ref_dir
            normal, plane_rmsd_val, plane_sv = fit_plane_normal(
                select_centered, reference_direction=orientation_ref,
            )

            ncl = float(plane_sv[1] / plane_sv[0] if plane_sv[0] > 1e-12 else 0.0)
            self._log_msg(
                f"[plane] fold={symmetry} rmsd={plane_rmsd_val:.4f}A "
                f"n_points={len(select_centered)} non_coll={ncl:.3f}"
            )
            if plane_rmsd_val > 2.0:
                self._log_msg(f"[plane] WARNING: plane RMSD {plane_rmsd_val:.2f}A exceeds 2.0A threshold")
            if ncl < 0.05:
                self._log_msg(
                    f"[plane] WARNING: near-collinear points (non_coll={ncl:.3f}); "
                    f"normal may be unreliable"
                )

            if global_axis is not None:
                dot_val = np.clip(np.abs(np.dot(normal, global_axis)), -1.0, 1.0)
                align_deg = float(np.rad2deg(np.arccos(dot_val)))
                self._log_msg(f"[plane] local_vs_global_axis_angle={align_deg:.2f}deg")
                if align_deg > 5.0:
                    self._log_msg(f"[plane] WARNING: alignment {align_deg:.2f}deg exceeds 5deg threshold")
                if align_deg > 10.0:
                    self._log_msg(
                        f"[plane] Overriding local normal with global axis "
                        f"(local alignment {align_deg:.1f}deg > 10deg)"
                    )
                    normal = global_axis

            if symmetry == 5 and len(select_centered) >= 5:
                penta = compute_pentagon_diagnostics(select_centered[:5], global_axis)
                ax_str = ""
                if 'alignment_angle_deg' in penta:
                    ax_str = f" vs_global_axis={penta['alignment_angle_deg']:.2f}deg"
                self._log_msg(
                    f"[pentagon] adj_dist={penta['adj_dist_mean']:.2f}±{penta['adj_dist_std']:.2f}A "
                    f"diag/adj={penta['diag_adj_ratio']:.3f} phi_dev={penta['phi_deviation']:.3f}{ax_str}"
                )

        z_axis = np.array([0.0, 0.0, 1.0])

        if rotation_check:
            R = rotation_matrix(normal, z_axis).T
            new_coords = np.asarray(coords_centered @ R, dtype=float)
            new_select_points = np.asarray(select_centered @ R, dtype=float)
        else:
            R = np.eye(3)
            new_coords = coords_centered.copy()
            new_select_points = select_centered.copy()

        self._emit_progress("align", 1.0)
        return AlignResult(
            aligned_coords=new_coords,
            rotation_matrix=R,
            select_points=new_select_points,
        )

    # ------------------------------------------------------------------
    # Stage 6 — Cleanup and format aligned coordinates
    # ------------------------------------------------------------------

    def format_aligned_pdb(
        self,
        pdb_info: PDBInfoResult,
        align_result: AlignResult,
    ) -> FormattedResult:
        """Center, min-zero, and PDB-format the aligned coordinates.

        Updates the chains DataFrame with formatted x/y/z/b/occ columns.
        """
        self._check_cancelled()
        self._emit_progress("format", 0.0)

        new_coords = align_result.aligned_coords.copy()

        # Center to COM
        new_coords = new_coords - np.mean(new_coords, axis=0)
        # Min-zero z
        new_coords = new_coords - np.array([0, 0, np.min(new_coords[:, 2])])
        # Min-zero x, y
        new_coords[:, 0] -= np.min(new_coords[:, 0])
        new_coords[:, 1] -= np.min(new_coords[:, 1])

        aligned_df = pd.DataFrame(new_coords, columns=["x", "y", "z"])

        chains = pdb_info.chains.copy()
        chains['x'] = aligned_df['x'].apply('{0:.3f}'.format)
        chains['y'] = aligned_df['y'].apply('{0:.3f}'.format)
        chains['z'] = aligned_df['z'].apply('{0:.3f}'.format)
        chains['b'] = chains['b'].astype('float').apply('{0:.2f}'.format)
        chains['occ'] = chains['occ'].astype('float').apply('{0:.2f}'.format)
        chains = chains.astype('string')
        chains = formatPDB.formatPDB(chains)

        self._emit_progress("format", 1.0)
        return FormattedResult(
            chains=chains,
            aligned_coords_df=aligned_df,
        )

    # ------------------------------------------------------------------
    # Stage 7 — Slice by weight
    # ------------------------------------------------------------------

    def slice_by_weight(
        self,
        formatted: FormattedResult,
        weight: float,
    ) -> SliceDataResult:
        """Apply the z-cutoff based on the slicing weight.

        ``cutoff = max(z) * weight``; atoms with z >= cutoff are kept.

        Returns
        -------
        SliceDataResult
        """
        self._check_cancelled()
        self._emit_progress("slice", 0.0)

        if not 0.0 <= weight <= 1.0:
            raise ValueError(f"Slicing weight must be between 0 and 1; got {weight}.")

        chains = formatted.chains
        z_vals = chains.z.astype(float)
        cutoff = np.max(z_vals) * weight

        # Record original chain residue counts for broken-chain detection
        sliced_chain_ids = np.unique(chains.chain[z_vals >= cutoff])
        original_lengths = []
        for x in sliced_chain_ids:
            original_lengths.append(
                len(np.unique(chains.resids[chains.chain == x]))
            )

        sliced_chains = chains[z_vals >= cutoff].copy()
        if len(sliced_chains) > 0:
            sliced_chains['z'] = (
                sliced_chains['z'].astype(float) - np.min(sliced_chains['z'].astype(float))
            )

        self._log_msg(f"[slice] cutoff_z={cutoff:.3f} weight={weight}")
        self._log_msg(f"[slice] atoms_before={len(chains)} atoms_after={len(sliced_chains)}")

        self._emit_progress("slice", 1.0)
        return SliceDataResult(
            sliced_chains=sliced_chains,
            cutoff=float(cutoff),
            original_lengths=original_lengths,
        )

    # ------------------------------------------------------------------
    # Stage 8 — Broken residue/chain cleanup
    # ------------------------------------------------------------------

    def cleanup_broken(
        self,
        sliced: SliceDataResult,
        amino_acid_checks: dict[str, int],
        orig_chains: list[str],
        broken_res_check: bool = True,
    ) -> CleanupResult:
        """Remove fragmented chains and broken boundary residues.

        Two-pass cleanup:
        1. Drop chains whose residue count is less than the original per-chain count.
        2. Drop residues near the slice edge (< 10 Å from bottom) that have
           fewer atoms than expected for their amino acid type.
        3. Re-run chain fragment check after residue removal.

        Returns
        -------
        CleanupResult
        """
        self._check_cancelled()
        self._emit_progress("cleanup", 0.0)

        sliced_chains = sliced.sliced_chains.copy()
        original_lengths = list(sliced.original_lengths)

        # --- Pass 1: Drop fragmented chains ---
        chain_values = np.unique(sliced_chains.chain)
        drop_list = []
        count = 0
        k = 0
        for i in tqdm(range(len(chain_values)), desc="Checking chain fragments"):
            self._check_cancelled()
            chain_val = chain_values[i]
            res_id = sliced_chains[sliced_chains.chain == chain_val].resids
            res_num = len(np.unique(res_id))
            if res_num < original_lengths[i]:
                drop_indices = sliced_chains[sliced_chains.chain == chain_val].index
                sliced_chains.drop(drop_indices, inplace=True)
                drop_list.append(i)
                count += 1
            k += 1
            if k == len(orig_chains):
                k = 0

        # Remove dropped entries from original_lengths
        i_off = 0
        for x in drop_list:
            del original_lengths[x - i_off]
            i_off += 1
        fragmented_chains_pass1 = count
        self._log_msg(f"Number of fragmented chains dropped: {count}")

        # --- Pass 2: Broken residue edge check ---
        fragmented_residues = 0
        if broken_res_check:
            edge_coords = sliced_chains[sliced_chains.z.astype(float) < 10.0]
            if len(edge_coords) > 0:
                chain_resid_values = np.vstack((
                    edge_coords.chain, edge_coords.resids,
                )).T
                check_values = np.vstack(tuple(set(map(tuple, chain_resid_values))))
                count = 0
                for i in tqdm(range(len(check_values)), desc="Checking edge residues"):
                    self._check_cancelled()
                    chain_val, resid_val = check_values[i]
                    resname = sliced_chains[
                        (sliced_chains.chain == chain_val)
                        & (sliced_chains.resids == resid_val)
                    ].resname
                    len_res = len(resname)
                    val_res = np.unique(resname)
                    if len(val_res) > 0:
                        resname_str = str(val_res[0]).strip()
                        len_check = amino_acid_checks.get(resname_str)
                        if len_check is not None and len_res != len_check:
                            drop_indices = sliced_chains[
                                (sliced_chains.chain == chain_val)
                                & (sliced_chains.resids == resid_val)
                            ].index
                            sliced_chains.drop(drop_indices, inplace=True)
                            count += 1
                fragmented_residues = count
                self._log_msg(f"Number of fragmented residues dropped: {count}")

                # Re-check for newly fragmented chains
                chain_values = np.unique(sliced_chains.chain)
                drop_list = []
                count = 0
                k = 0
                for i in tqdm(range(len(chain_values)), desc="Re-checking chain fragments"):
                    self._check_cancelled()
                    chain_val = chain_values[i]
                    res_id = sliced_chains[sliced_chains.chain == chain_val].resids
                    res_num = len(np.unique(res_id))
                    if res_num < original_lengths[i]:
                        drop_indices = sliced_chains[sliced_chains.chain == chain_val].index
                        sliced_chains.drop(drop_indices, inplace=True)
                        drop_list.append(i)
                        count += 1
                    k += 1
                    if k == len(orig_chains):
                        k = 0
                self._log_msg(f"Number of fragmented chains dropped: {count}")
                i_off = 0
                for x in drop_list:
                    del original_lengths[x - i_off]
                    i_off += 1
                fragmented_chains_pass2 = count
            else:
                self._log_msg("No fragmented residues. Nice!")
                fragmented_chains_pass2 = 0
        else:
            fragmented_chains_pass2 = 0

        self._emit_progress("cleanup", 1.0)
        return CleanupResult(
            cleaned_chains=sliced_chains,
            fragmented_chains_dropped=fragmented_chains_pass1 + fragmented_chains_pass2,
            fragmented_residues_dropped=fragmented_residues,
            original_lengths=original_lengths,
        )

    # ------------------------------------------------------------------
    # Stage 9 — Rename chains and write output
    # ------------------------------------------------------------------

    def finalize_and_write(
        self,
        cleaned: CleanupResult,
        output_path: str,
    ) -> SliceOutput:
        """Rename chains, re-index atoms, format PDB, and write to disk.

        Parameters
        ----------
        cleaned : CleanupResult
        output_path : str
            Absolute or relative path for the output PDB file.

        Returns
        -------
        SliceOutput
        """
        self._check_cancelled()
        self._emit_progress("write", 0.0)

        sliced_chains = cleaned.cleaned_chains.copy()

        if len(sliced_chains) == 0:
            out_dir = os.path.dirname(output_path)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            with open(output_path, "w", encoding="utf-8"):
                pass
            self._log_msg(f"Saved empty sliced capsid to {output_path}")
            self._emit_progress("write", 1.0)
            return SliceOutput(
                output_pdb_path=output_path,
                output_atoms=0,
                output_residues=0,
                output_chains=0,
                log_lines=list(self._log),
            )

        # Rename chains
        chain_values = np.unique(sliced_chains.chain)
        chain_dict = createDictionary.createChainDictionary(chain_values)
        sliced_chains.chain = sliced_chains.chain.replace(chain_dict)

        chain_list = list(np.unique(sliced_chains.chain))
        self._log_msg(f"Output chains: {chain_list}")
        if len(chain_list) > 60:
            self._log_msg("Heads up! More than 60 unique chains. Some IDs may not be recognized.")

        # Re-base z
        sliced_chains["z"] = (
            sliced_chains.z.astype(float) - np.min(sliced_chains.z.astype(float))
        )

        # Re-index atoms
        sliced_chains["idx"] = range(1, len(sliced_chains) + 1)

        # Count input (before cleanup happens upstream; we track from the full chains)
        input_atoms = len(cleaned.cleaned_chains)  # This is post-cleanup — we need pre-cleanup
        # We won't have pre-cleanup here, but we can compute from the sliced chains before drop

        # Final formatting and save
        unformatted = copy.deepcopy(sliced_chains)
        # viewer_chain_index is internal metadata and must never be serialized
        # into the fixed-column scientific PDB output.
        pdb_chains = sliced_chains.drop(
            columns=["viewer_chain_index"],
            errors="ignore",
        )
        formatted = formatPDB.formatPDB(pdb_chains)

        # Ensure output directory exists
        out_dir = os.path.dirname(output_path)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        try:
            savetxt(output_path, formatted, fmt='%s', delimiter='')
        except (OSError, ValueError) as e:
            raise RuntimeError(f"Couldn't save {output_path}: {e}") from e

        self._log_msg(f"Saved sliced capsid to {output_path}")

        # Compute statistics
        output_atoms = len(formatted)
        output_residues = int(
            unformatted[["chain", "resids"]].drop_duplicates().shape[0]
        )
        output_chains = len(chain_list)

        self._emit_progress("write", 1.0)
        return SliceOutput(
            output_pdb_path=output_path,
            input_atoms=0,  # populated by run_full
            input_residues=0,
            input_chains=0,
            output_atoms=output_atoms,
            output_residues=output_residues,
            output_chains=output_chains,
            retained_percentage=0.0,  # populated by run_full
            axis_dir=[0.0, 0.0, 0.0],
            diagnostics={},
            log_lines=list(self._log),
        )

    # ------------------------------------------------------------------
    # Composite: run_full
    # ------------------------------------------------------------------

    def run_full(self, request: SliceRequest) -> SliceOutput:
        """Run the complete pipeline end-to-end.

        Parameters
        ----------
        request : SliceRequest
            All parameters for the slice operation.

        Returns
        -------
        SliceOutput
        """
        # Stage 1: Load
        load_result = self.load_structure(request.input_path)

        # Stage 2: Detect axis
        axis_result = self.detect_axis(
            symmetry=request.symmetry,
            auto=request.auto,
            axis_index=request.axis_index,
            roi_selection=request.roi_selection,
            roi_frame=request.roi_frame,
            ref_indices=request.ref_indices,
            legacy_plane_heuristic=request.legacy_plane_heuristic,
        )

        # Stage 3: Select plane points
        plane_result = self.select_plane_points(
            axis_result=axis_result,
            symmetry=request.symmetry,
            auto_mode=request.auto,
            legacy_plane_heuristic=request.legacy_plane_heuristic,
            plot=request.plot,
        )

        # Stage 4: Extract PDB info
        pdb_info = self.extract_pdb_info()

        # Stage 5: Align
        align_result = self.align_coordinates(
            pdb_info=pdb_info,
            plane_result=plane_result,
            symmetry=request.symmetry,
            auto_mode=request.auto,
            legacy_plane_heuristic=request.legacy_plane_heuristic,
            plot=request.plot,
        )

        # Stage 6: Format aligned PDB
        formatted = self.format_aligned_pdb(
            pdb_info=pdb_info,
            align_result=align_result,
        )

        # Record pre-slice counts
        input_atoms = len(formatted.chains)
        input_residues = load_result.total_residues
        input_chains = len(np.unique(pdb_info.chains.chain))

        # Stage 7: Slice
        sliced = self.slice_by_weight(formatted=formatted, weight=request.weight)

        # Stage 8: Cleanup broken
        cleaned = self.cleanup_broken(
            sliced=sliced,
            amino_acid_checks=pdb_info.amino_acid_checks,
            orig_chains=load_result.orig_chains,
        )

        # Stage 9: Finalize and write
        output_path = request.workspace_path or ""
        if output_path:
            output_path = os.path.join(
                output_path,
                f"{os.path.splitext(os.path.basename(request.input_path))[0]}"
                f"-sliced-sym{request.symmetry}-w{request.weight}.pdb",
            )
        else:
            output_path = (
                f"output/{os.path.splitext(os.path.basename(request.input_path))[0]}"
                f"-sliced-sym{request.symmetry}-w{request.weight}.pdb"
            )

        result = self.finalize_and_write(cleaned=cleaned, output_path=output_path)

        # Preserve the original flattened chain IDs for the viewer so it can
        # map retained atoms back to M{frame}_{chain}. Match the final PDB's
        # z rebasing without applying its scientific-output chain renaming.
        viewer_chains = cleaned.cleaned_chains.copy()
        if len(viewer_chains) > 0:
            viewer_chains["z"] = (
                viewer_chains.z.astype(float)
                - np.min(viewer_chains.z.astype(float))
            )
        result.viewer_chains = viewer_chains

        # Fill in pre-slice counts
        result.input_atoms = input_atoms
        result.input_residues = input_residues
        result.input_chains = input_chains
        result.retained_percentage = (
            (result.output_atoms / input_atoms * 100.0) if input_atoms > 0 else 0.0
        )
        result.axis_dir = axis_result.axis_dir
        result.diagnostics = {
            "plane_rmsd": plane_result.plane_rmsd,
            "non_collinearity": plane_result.non_collinearity,
            "local_vs_global_angle_deg": plane_result.local_vs_global_angle_deg,
            "fragmented_chains_dropped": cleaned.fragmented_chains_dropped,
            "fragmented_residues_dropped": cleaned.fragmented_residues_dropped,
            "cutoff": sliced.cutoff,
        }

        return result

    # ------------------------------------------------------------------
    # Composite: run_cli (backward-compatible with sliceCapsid.py CLI)
    # ------------------------------------------------------------------

    @staticmethod
    def run_cli(args: Any) -> int:
        """Run the pipeline with argparse Namespace, preserving legacy CLI behavior.

        Parameters
        ----------
        args : argparse.Namespace
            Parsed CLI arguments (matching sliceCapsid.py's argument definitions).

        Returns
        -------
        int
            Exit code (0 on success, 1 on error).
        """
        symmetry = args.symmetry
        if symmetry not in (2, 3, 5):
            raise ValueError(f"Unsupported symmetry {symmetry}; expected one of 2, 3, or 5.")

        auto_mode = args.auto
        roi_selection = args.roi_selection
        roi_frame = args.roi_frame
        weight = args.weight
        axis_index = args.axis_index
        legacy_plane = args.legacy_plane_heuristic
        output_title = args.pdb
        ref_indices = args.refindex

        # Input paths (preserve legacy relative path behavior)
        pdb_path = f"input/{output_title}.pdb"

        # Handle --list-axes without --auto (print and exit)
        if args.list_axes and not auto_mode:
            u = md.Universe(pdb_path, pdb_path)
            axis_rows = findSymmetryAxes.list_axes(
                u, symmetry,
                roi_selection=roi_selection,
                roi_frame=roi_frame,
            )
            _print_axis_table(axis_rows, has_roi=bool(roi_selection))
            return 0

        if roi_selection and not auto_mode:
            print("[roi] --roi-selection is only used with --auto or --list-axes.")

        # Build request and run
        request = SliceRequest(
            input_path=pdb_path,
            workspace_path="output",
            symmetry=symmetry,
            weight=weight,
            auto=auto_mode,
            axis_index=axis_index,
            roi_selection=roi_selection,
            roi_frame=roi_frame,
            ref_indices=ref_indices,
            legacy_plane_heuristic=legacy_plane,
            plot=args.plot,
        )

        pipeline = CapsidPipeline()
        result = pipeline.run_full(request)

        # Print log messages (preserve legacy stdout behavior)
        for msg in result.log_lines:
            print(msg)
        print(f"Retained: {result.retained_percentage:.1f}% "
              f"({result.output_atoms}/{result.input_atoms} atoms)")

        return 0


# ---------------------------------------------------------------------------
# Legacy helper (preserved for CLI backward compat)
# ---------------------------------------------------------------------------


def _print_axis_table(rows: list, has_roi: bool = False) -> None:
    """Print axis candidates in the legacy table format."""
    if len(rows) == 0:
        print("[axes] No axes found.")
        return

    if has_roi:
        print("[axes] rank axis-index roi-score roi-angle-deg roi-line-dist symmetry-score axis")
        for rank, row in enumerate(rows):
            axis = row["axis"]
            print(
                "[axes] "
                f"{rank:>4} {row['axis_index']:>10} "
                f"{row['roi_score']:.5f} {row['roi_angle_deg']:.3f} "
                f"{row['roi_line_distance']:.3f} {row['score']:.5f} "
                f"[{axis[0]:.5f}, {axis[1]:.5f}, {axis[2]:.5f}]"
            )
    else:
        print("[axes] rank axis-index symmetry-score ref-index ref-frame chain axis")
        for rank, row in enumerate(rows):
            axis = row["axis"]
            print(
                "[axes] "
                f"{rank:>4} {row['axis_index']:>10} {row['score']:.5f} "
                f"{row['ref_index']:>9} {row['ref_frame']:>9} "
                f"{row['chain_id']} "
                f"[{axis[0]:.5f}, {axis[1]:.5f}, {axis[2]:.5f}]"
            )
