#!/usr/bin/env python3
"""
Automated Icosahedral Symmetry Axis Detection for CapSACIN.

Given an icosahedral virus capsid structure (MDAnalysis Universe), this module
automatically detects the 5-fold, 3-fold, and 2-fold symmetry axes using the
distribution of protein chain centers of mass, then selects a reference atom
near a chosen axis — eliminating the need for users to manually provide
reference atom indices via VMD.

Method (Direct Spherical Search):
    1. Compute one COM per MODEL/asymmetric-unit copy and the capsid center.
    2. Sample directions on the unit sphere with a Fibonacci grid.
    3. Score each direction by how well rotations by 360/fold map monomer
       directions onto other monomer directions.
    4. Deduplicate local minima to obtain candidate 5-fold, 3-fold, or 2-fold
       axes, ranked by symmetry score.
    5. Select the requested axis and choose representative reference atoms.

Reference
---------
The axis formulas follow the standard icosahedral coordinates in a basis where
three orthogonal 2-fold axes are aligned with the coordinate axes.
"""

import numpy as np
import MDAnalysis as md
from MDAnalysis.analysis import contacts
from capsacin.definePlane import fit_plane_normal, select_5fold_ring


# ---------------------------------------------------------------------------
# Golden ratio
# ---------------------------------------------------------------------------
PHI = (1.0 + np.sqrt(5.0)) / 2.0  # ≈ 1.6180339887
PHI_INV = 1.0 / PHI               # ≈ 0.6180339887


# ---------------------------------------------------------------------------
# Step 1: Capsid centre and chain centres of mass
# ---------------------------------------------------------------------------

def compute_capsid_center(universe):
    """
    Compute the centre of mass of all protein atoms across **all frames**.

    For a multi-MODEL capsid PDB, this averages over all 60 icosahedral copies
    to obtain the capsid-wide centre.

    Parameters
    ----------
    universe : MDAnalysis.Universe

    Returns
    -------
    center : np.ndarray of shape (3,)
    """
    all_positions = []
    for ts in universe.trajectory:
        protein = universe.select_atoms("protein")
        all_positions.append(protein.positions)
    stacked = np.vstack(all_positions)
    return stacked.mean(axis=0)


def compute_monomer_coms(universe):
    """
    Compute the centre of mass for each monomer (asymmetric unit).

    CapSACIN input PDBs store the 60 icosahedral copies as separate MODELs
    (NMR-style multi-model format), where each frame is one copy of the
    asymmetric unit.  This function iterates over trajectory frames and
    returns one COM per frame.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        A multi-model PDB with one asymmetric unit per MODEL.

    Returns
    -------
    monomer_coms : np.ndarray of shape (N_frames, 3)
        COM of each monomer (one per trajectory frame).
    """
    coms = []
    for ts in universe.trajectory:
        protein = universe.select_atoms("protein")
        coms.append(protein.positions.mean(axis=0))
    return np.array(coms)


# ---------------------------------------------------------------------------
# Step 2: Direct spherical search for rotational symmetry axes
# ---------------------------------------------------------------------------
# Because the covariance matrix eigenvectors are arbitrary for a nearly-
# isotropic icosahedral distribution, we use a direct grid search over the
# sphere to find directions that maximise n-fold rotational symmetry.

def _fibonacci_sphere(n_samples):
    """
    Generate *n_samples* approximately uniform points on the unit sphere
    using the Fibonacci lattice method.

    Returns
    -------
    directions : np.ndarray of shape (n_samples, 3)
    """
    if n_samples < 2:
        raise ValueError(f"n_samples must be >= 2, got {n_samples}.")
    golden = (1.0 + np.sqrt(5.0)) / 2.0
    dirs = np.empty((n_samples, 3))
    for i in range(n_samples):
        y = 1.0 - (i / float(n_samples - 1)) * 2.0  # [1, -1]
        radius = np.sqrt(1.0 - y * y)
        theta = 2.0 * np.pi * i / golden
        dirs[i, 0] = np.cos(theta) * radius
        dirs[i, 1] = y
        dirs[i, 2] = np.sin(theta) * radius
    return dirs


def _rotational_symmetry_score(axis, monomer_dirs, fold):
    """
    Score how well *axis* acts as an *fold*-fold rotational symmetry axis.

    Vectorised over monomers for speed.

    Parameters
    ----------
    axis : np.ndarray of shape (3,)
        Candidate axis direction (need not be unit).
    monomer_dirs : np.ndarray of shape (N, 3)
        Unit vectors from capsid centre to each monomer COM.
    fold : int
        Rotational symmetry order (2, 3, or 5).

    Returns
    -------
    score : float
        Mean angular error in radians.  0 = perfect *fold*-fold symmetry.
    """
    axis = _normalize(axis)
    N = len(monomer_dirs)
    total_err = 0.0

    for k in range(1, fold):
        angle = k * 2.0 * np.pi / fold
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        # Rodrigues rotation of ALL monomer directions at once
        dot_products = monomer_dirs @ axis  # (N,)
        cross_products = np.cross(axis, monomer_dirs)  # (N, 3)
        rotated = (
            cos_a * monomer_dirs
            + sin_a * cross_products
            + (1.0 - cos_a) * np.outer(dot_products, axis)
        )
        # Normalise all rotated vectors
        norms = np.linalg.norm(rotated, axis=1, keepdims=True)
        rotated = rotated / norms

        # For each rotated monomer, find nearest actual monomer
        # cos_dists is (N, N): cos_dists[i,j] = |rotated[i] · monomer[j]|
        cos_dists = np.abs(rotated @ monomer_dirs.T)  # (N, N)
        cos_dists = np.clip(cos_dists, -1.0, 1.0)
        min_angles = np.min(np.arccos(cos_dists), axis=1)  # (N,)
        total_err += np.sum(min_angles)

    return total_err / (N * (fold - 1))


def _find_local_minima(scores, directions, min_angle=np.deg2rad(5.0)):
    """
    Find local minima in the score landscape on the sphere.

    A direction is a local minimum if its score is lower than all neighbours
    within *min_angle*.

    Parameters
    ----------
    scores : np.ndarray of shape (N,)
    directions : np.ndarray of shape (N, 3)
    min_angle : float
        Angular radius for neighbourhood check (radians).

    Returns
    -------
    minima : list of (score, direction) tuples, sorted by ascending score.
    """
    cos_thresh = np.cos(min_angle)
    minima = []
    for i in range(len(scores)):
        # Neighbours: all points within min_angle
        cos_dists = np.abs(directions @ directions[i])
        neighbours = np.where(cos_dists > cos_thresh)[0]
        if all(scores[i] <= scores[j] for j in neighbours):
            minima.append((scores[i], directions[i].copy()))
    minima.sort(key=lambda x: x[0])
    return minima


def _deduplicate_minima(minima, min_separation=np.deg2rad(10.0)):
    """
    Remove duplicate minima that are closer than *min_separation*.

    Parameters
    ----------
    minima : list of (score, direction)
    min_separation : float
        Minimum angular separation between distinct axes (radians).

    Returns
    -------
    unique : list of (score, direction)
    """
    cos_sep = np.cos(min_separation)
    unique = []
    for score, direction in minima:
        is_dup = False
        for _, existing in unique:
            if abs(np.dot(direction, existing)) > cos_sep:
                is_dup = True
                break
        if not is_dup:
            unique.append((score, direction))
    return unique


def find_symmetry_axes(monomer_coms, capsid_center, fold,
                       n_grid=5000, n_expected=None):
    """
    Find *fold*-fold rotational symmetry axes via direct spherical search.

    Parameters
    ----------
    monomer_coms : np.ndarray of shape (N, 3)
        Monomer COM positions.
    capsid_center : np.ndarray of shape (3,)
        Capsid centre.
    fold : int
        Rotational symmetry order: 2, 3, or 5.
    n_grid : int
        Number of sample points on the sphere.
    n_expected : int or None
        Expected number of axes (6 for 5-fold, 10 for 3-fold, 15 for
        2-fold).  If None, all detected axes are returned.

    Returns
    -------
    axes : np.ndarray of shape (K, 3)
        Detected axis directions (unit vectors), sorted by ascending
        symmetry score (best first).
    """
    monomer_dirs = monomer_coms - capsid_center
    monomer_dirs = np.array([_normalize(d) for d in monomer_dirs])

    # Generate sample directions on the sphere
    directions = _fibonacci_sphere(n_grid)

    # Score each direction
    scores = np.array([
        _rotational_symmetry_score(d, monomer_dirs, fold)
        for d in directions
    ])

    # Find local minima
    minima = _find_local_minima(scores, directions)

    # Deduplicate
    unique = _deduplicate_minima(minima)

    # Take the n_expected best
    if n_expected is not None:
        unique = unique[:n_expected]

    axes = np.array([d for _, d in unique])
    return axes


def build_icosahedral_axes(monomer_coms, capsid_center, n_grid=2000):
    """
    Find all icosahedral symmetry axes via direct spherical search.

    Parameters
    ----------
    monomer_coms : np.ndarray of shape (N, 3)
    capsid_center : np.ndarray of shape (3,)
    n_grid : int
        Number of sample points on the sphere.

    Returns
    -------
    axes : dict
        {'5fold': (K5, 3), '3fold': (K3, 3), '2fold': (K2, 3)}
    """
    fivefold = find_symmetry_axes(
        monomer_coms, capsid_center, fold=5,
        n_grid=n_grid, n_expected=6,
    )
    threefold = find_symmetry_axes(
        monomer_coms, capsid_center, fold=3,
        n_grid=n_grid, n_expected=10,
    )
    twofold = find_symmetry_axes(
        monomer_coms, capsid_center, fold=2,
        n_grid=n_grid, n_expected=15,
    )
    return {"5fold": fivefold, "3fold": threefold, "2fold": twofold}


def _axis_key(symmetry_type):
    if symmetry_type not in (2, 3, 5):
        raise ValueError(f"Invalid symmetry_type: {symmetry_type}")
    return {2: "2fold", 3: "3fold", 5: "5fold"}[symmetry_type]


def _compute_axis_context(universe, n_grid=2000):
    capsid_center = compute_capsid_center(universe)
    monomer_coms = compute_monomer_coms(universe)

    if len(monomer_coms) < 3:
        raise ValueError(
            f"Found only {len(monomer_coms)} monomers (MODELs). "
            "Need at least 3 for icosahedral symmetry detection."
        )

    monomer_dirs = monomer_coms - capsid_center
    monomer_dirs = np.array([_normalize(d) for d in monomer_dirs])
    all_axes = build_icosahedral_axes(monomer_coms, capsid_center, n_grid=n_grid)

    return capsid_center, monomer_coms, monomer_dirs, all_axes


def compute_roi_center(universe, roi_selection, roi_frame=0):
    """
    Compute the geometric center of a user-defined ROI in one capsid copy.

    ``roi_frame`` defaults to 0 to match the legacy ``--refindex`` behavior:
    the user describes the ROI on one physical copy, and symmetry operations
    find its related copies elsewhere in the multi-MODEL capsid.
    """
    if not roi_selection:
        raise ValueError("roi_selection must be provided.")
    if roi_frame < 0 or roi_frame >= universe.trajectory.n_frames:
        raise ValueError(
            f"roi_frame={roi_frame} but the trajectory has "
            f"{universe.trajectory.n_frames} frames."
        )

    universe.trajectory[roi_frame]
    atoms = universe.select_atoms(roi_selection)
    if len(atoms) == 0:
        raise ValueError(
            f"ROI selection matched no atoms on frame {roi_frame}: {roi_selection!r}"
        )
    center = atoms.positions.mean(axis=0)
    universe.trajectory[0]
    return center


def _rank_axes_for_roi_from_context(
    candidates,
    monomer_coms,
    monomer_dirs,
    capsid_center,
    roi_center,
    symmetry_type,
):
    roi_vector = roi_center - capsid_center
    roi_distance = np.linalg.norm(roi_vector)
    if roi_distance < 1e-12:
        raise ValueError("ROI center is too close to the capsid center.")

    roi_dir = roi_vector / roi_distance
    capsid_radius = np.mean(np.linalg.norm(monomer_coms - capsid_center, axis=1))
    if capsid_radius < 1e-12:
        capsid_radius = 1.0

    results = []
    for axis_index, axis in enumerate(candidates):
        axis = _normalize(axis)
        if np.dot(axis, roi_dir) < 0:
            axis = -axis

        axis_dot = np.clip(np.dot(axis, roi_dir), -1.0, 1.0)
        roi_angle = np.arccos(axis_dot)
        roi_line_distance = np.linalg.norm(np.cross(roi_vector, axis))
        roi_distance_score = roi_line_distance / capsid_radius
        symmetry_score = _rotational_symmetry_score(
            axis, monomer_dirs, symmetry_type
        )

        # ROI proximity is the main selector; symmetry score breaks ties among
        # equivalent axes and keeps noisy structures from picking poor minima.
        roi_score = roi_angle + roi_distance_score + 0.25 * symmetry_score

        results.append({
            "axis_index": axis_index,
            "axis": axis,
            "roi_score": roi_score,
            "roi_angle_deg": np.rad2deg(roi_angle),
            "roi_line_distance": roi_line_distance,
            "roi_distance_score": roi_distance_score,
            "symmetry_score": symmetry_score,
        })

    results.sort(key=lambda x: x["roi_score"])
    return results


def rank_axes_for_roi(universe, symmetry_type, roi_selection, roi_frame=0,
                      n_grid=2000):
    """
    Rank candidate symmetry axes by proximity to a user-defined ROI.

    Returns a list of dictionaries sorted by increasing ``roi_score``.  The
    returned axis vectors are oriented so their positive direction points
    toward the ROI, which is important because CapSACIN keeps the high-z side
    after alignment.
    """
    capsid_center, monomer_coms, monomer_dirs, all_axes = _compute_axis_context(
        universe, n_grid=n_grid
    )
    roi_center = compute_roi_center(universe, roi_selection, roi_frame=roi_frame)
    key = _axis_key(symmetry_type)

    return _rank_axes_for_roi_from_context(
        all_axes[key],
        monomer_coms,
        monomer_dirs,
        capsid_center,
        roi_center,
        symmetry_type,
    )


# Backward-compatible alias used by auto_detect_reference
def find_principal_axes(monomer_coms, capsid_center):
    """
    Deprecated — kept for API compatibility.
    Use :func:`find_symmetry_axes` with ``fold=2`` instead.
    """
    return find_symmetry_axes(
        monomer_coms, capsid_center, fold=2, n_grid=5000, n_expected=3
    ).T  # (3,3) format


# ---------------------------------------------------------------------------
# Utility: vector normalisation
# ---------------------------------------------------------------------------

def _normalize(v):
    """Return a unit vector in the direction of *v*."""
    n = np.linalg.norm(v)
    if n < 1e-12:
        return v
    return v / n


# ---------------------------------------------------------------------------
# Step 5: Find a reference atom near the chosen axis
# ---------------------------------------------------------------------------

def find_reference_atom(universe, axis, capsid_center, symmetry_type):
    """
    Find the best reference atom for the given symmetry axis.

    Iterates over candidate atoms in the monomer nearest *axis*, computes the
    plane normal from its 60 symmetry-related copies using the same
    nearest-neighbour distance-ranking heuristic as ``sliceCapsid.py``, and
    returns the atom whose normal best aligns with *axis*.

    Parameters
    ----------
    universe : MDAnalysis.Universe
    axis : np.ndarray of shape (3,)
        Symmetry axis direction (unit vector).
    capsid_center : np.ndarray of shape (3,)
        Capsid centre.
    symmetry_type : int
        2, 3, or 5.

    Returns
    -------
    ref_index : int
        0-based atom index suitable for use as ``--refindex``.
    ref_frame : int
        The trajectory frame (0-based) of the monomer closest to the axis.
    """
    monomer_coms = compute_monomer_coms(universe)
    monomer_dirs = monomer_coms - capsid_center
    monomer_dirs = np.array([_normalize(d) for d in monomer_dirs])

    # Find the monomer closest to the axis
    cos_angles = np.abs(monomer_dirs @ axis)
    cos_angles = np.clip(cos_angles, -1.0, 1.0)
    angles = np.arccos(cos_angles)
    best_frame = int(np.argmin(angles))

    # Collect all atom positions from all frames for the distance-based search
    # (same approach as sliceCapsid.py)
    all_frame_coords = []
    all_frame_atoms = []
    for ts in universe.trajectory:
        protein = universe.select_atoms("protein")
        all_frame_coords.append(protein.positions.copy())
        # store atom metadata for matching
        all_frame_atoms.append({
            "resnames": protein.resnames.copy(),
            "resids": protein.resids.copy(),
            "names": protein.names.copy(),
            "indices": protein.indices.copy(),
        })

    # Determine the nearest-neighbour ranks for this symmetry type
    # (matching sliceCapsid.py logic)
    if symmetry_type == 3:
        k1, k2 = 1, 2
    elif symmetry_type == 5:
        k1, k2 = 1, 3
    elif symmetry_type == 2:
        k1, k2 = 0, 1
    else:
        raise ValueError(f"Invalid symmetry_type: {symmetry_type}")

    # Candidate atoms: sample Cα atoms plus surface-exposed atoms from the
    # best-frame monomer.  Testing every atom is cheap enough (~4000 atoms).
    universe.trajectory[best_frame]
    monomer_atoms = universe.select_atoms("protein")
    n_atoms = len(monomer_atoms)

    # Build all-copies coordinate arrays for fast distance computation
    # For each atom type (resname+resid+name), we need to find all copies.
    # This is expensive to do per atom, so we pre-build.
    #
    # Strategy: for each candidate atom *i* in the best frame, find its
    # 60 copies across all frames, compute distances from copy in best_frame,
    # pick k1/k2 nearest, compute plane normal, score vs axis.
    #
    # We sample Cα atoms and a few other atoms for good coverage.
    # Cα atoms trace the protein backbone and are well-distributed.
    ca_indices = [
        i for i in range(n_atoms) if monomer_atoms.names[i] == "CA"
    ]
    # Also sample every 10th non-CA atom for surface coverage
    other_indices = [
        i for i in range(n_atoms)
        if monomer_atoms.names[i] != "CA" and i % 10 == 0
    ]
    candidate_local_indices = ca_indices + other_indices
    # Limit to ~600 candidates max
    if len(candidate_local_indices) > 600:
        candidate_local_indices = candidate_local_indices[:600]

    best_score = -np.inf
    best_idx = monomer_atoms.indices[0]
    best_normal = None

    for local_i in candidate_local_indices:
        ref_idx_in_frame = monomer_atoms.indices[local_i]
        ref_resname = monomer_atoms.resnames[local_i]
        ref_resid = monomer_atoms.resids[local_i]
        ref_name = monomer_atoms.names[local_i]

        # Find all copies across all frames
        copies = []
        for frame_idx in range(len(all_frame_coords)):
            af = all_frame_atoms[frame_idx]
            mask = (
                (af["resnames"] == ref_resname)
                & (af["resids"] == ref_resid)
                & (af["names"] == ref_name)
            )
            if np.any(mask):
                copies.append(all_frame_coords[frame_idx][mask][0])

        min_copies = 5 if symmetry_type == 5 else max(k1, k2) + 1
        if len(copies) < min_copies:
            continue

        copies = np.array(copies)

        # The copy in best_frame is pointA
        # Find it among copies (the one with matching frame and index)
        pointA = all_frame_coords[best_frame][
            (all_frame_atoms[best_frame]["resnames"] == ref_resname)
            & (all_frame_atoms[best_frame]["resids"] == ref_resid)
            & (all_frame_atoms[best_frame]["names"] == ref_name)
        ][0]

        # Distances from pointA to all copies
        dists = np.linalg.norm(copies - pointA, axis=1)
        dist_sort = np.argsort(dists)

        # Compute plane normal via SVD-based fit
        if symmetry_type == 5:
            # Use all 5 ring copies for a better plane fit
            ring_points = select_5fold_ring(copies, pointA)
            normal, _rmsd, sv = fit_plane_normal(ring_points,
                                                  reference_direction=axis)
            # For backward-compat pointB/pointC references
            pointB = ring_points[1].copy()
            pointC = ring_points[3].copy()
        else:
            # 3-fold: use 3 points directly; 2-fold: keep z-perturbation
            # (z-perturbation is adequate for candidate scoring)
            pointB = copies[dist_sort[k1]].copy()
            pointC = copies[dist_sort[k2]].copy()
            if symmetry_type == 2:
                pointB[2] += 0.1
                pointC[2] -= 0.1
            normal, _rmsd, sv = fit_plane_normal(
                np.array([pointA, pointB, pointC]),
                reference_direction=axis,
            )

        # Score: how well does this normal align with the theoretical axis?
        alignment = abs(np.dot(normal, axis))
        # Normalized non-collinearity: s[1]/s[0] ∈ [0,1]
        #   → near 1: well-spread points, normal is reliable
        #   → near 0: nearly collinear, normal is noise
        ncl = sv[1] / sv[0] if sv[0] > 1e-12 else 0.0
        score = alignment + 0.05 * min(ncl, 1.0)

        if score > best_score:
            best_score = score
            best_idx = ref_idx_in_frame
            best_normal = normal

    # Restore frame 0 for compatibility
    universe.trajectory[0]

    return best_idx, best_frame


# ---------------------------------------------------------------------------
# Top-level API
# ---------------------------------------------------------------------------

def _find_reference_frames(monomer_dirs, best_axis, symmetry_type,
                           prefer_positive_axis=False):
    # Pick monomer A.  ROI-aware axes are oriented toward the ROI, so their
    # positive direction matters.  Without an ROI, preserve the legacy behavior
    # where the two directions of one physical axis are equivalent.
    cos_angles = monomer_dirs @ best_axis
    cos_angles = np.clip(cos_angles, -1.0, 1.0)
    if prefer_positive_axis:
        frame_a = int(np.argmax(cos_angles))
    else:
        frame_a = int(np.argmax(np.abs(cos_angles)))

    angle_step = 2.0 * np.pi / symmetry_type

    rotated_180 = _normalize(
        np.cos(np.pi) * monomer_dirs[frame_a]
        + np.sin(np.pi) * np.cross(best_axis, monomer_dirs[frame_a])
        + (1.0 - np.cos(np.pi)) * np.dot(best_axis, monomer_dirs[frame_a]) * best_axis
    )
    cos_dists_180 = monomer_dirs @ rotated_180
    if not prefer_positive_axis:
        cos_dists_180 = np.abs(cos_dists_180)
    cos_dists_180 = np.clip(cos_dists_180, -1.0, 1.0)

    rotated_dir = _normalize(
        np.cos(angle_step) * monomer_dirs[frame_a]
        + np.sin(angle_step) * np.cross(best_axis, monomer_dirs[frame_a])
        + (1.0 - np.cos(angle_step)) * np.dot(best_axis, monomer_dirs[frame_a]) * best_axis
    )
    cos_dists = monomer_dirs @ rotated_dir
    if not prefer_positive_axis:
        cos_dists = np.abs(cos_dists)
    cos_dists = np.clip(cos_dists, -1.0, 1.0)
    frame_b = int(np.argmax(cos_dists))

    if symmetry_type == 2:
        frame_b = frame_a
        frame_c = int(np.argmax(cos_dists_180))
    else:
        rotated_dir2 = _normalize(
            np.cos(2 * angle_step) * monomer_dirs[frame_a]
            + np.sin(2 * angle_step) * np.cross(best_axis, monomer_dirs[frame_a])
            + (1.0 - np.cos(2 * angle_step)) * np.dot(best_axis, monomer_dirs[frame_a]) * best_axis
        )
        cos_dists2 = monomer_dirs @ rotated_dir2
        if not prefer_positive_axis:
            cos_dists2 = np.abs(cos_dists2)
        cos_dists2 = np.clip(cos_dists2, -1.0, 1.0)
        frame_c = int(np.argmax(cos_dists2))

    return [frame_a, frame_b, frame_c]


def find_symmetry_related_frames(universe, axis, symmetry_type, seed_frame=None):
    """Return the MODEL frames forming one local n-fold feature.

    The multi-MODEL capsid inputs store one symmetry copy per trajectory
    frame. When ``seed_frame`` is omitted, the monomer nearest the selected
    physical axis is chosen automatically. This is the correct behavior for
    globally ranked axes. When ROI-aware ranking is active, the caller can
    provide the user's ROI frame as the seed. Rotating that monomer direction
    by ``360 / symmetry_type`` degrees identifies the complete local feature.
    """
    if symmetry_type not in (2, 3, 5):
        raise ValueError(
            f"symmetry_type must be 2, 3, or 5; got {symmetry_type}."
        )

    n_frames = universe.trajectory.n_frames
    if seed_frame is not None and (seed_frame < 0 or seed_frame >= n_frames):
        raise ValueError(
            f"seed_frame={seed_frame} but the trajectory has {n_frames} frames."
        )

    center = compute_capsid_center(universe)
    monomer_coms = compute_monomer_coms(universe)
    centered_coms = monomer_coms - center
    norms = np.linalg.norm(centered_coms, axis=1, keepdims=True)
    monomer_dirs = centered_coms / np.maximum(norms, 1e-12)
    axis = _normalize(np.asarray(axis, dtype=float))
    if seed_frame is None:
        # The alignment maps this signed axis to +z and slicing retains the
        # high-z cap. Pick the local feature at that same (+axis) end so its
        # ROI remains visible on the sliced/overlay structure after Run.
        seed_frame = int(np.argmax(monomer_dirs @ axis))
    seed_dir = monomer_dirs[seed_frame]

    frames = []
    available = np.ones(n_frames, dtype=bool)
    angle_step = 2.0 * np.pi / symmetry_type

    for step in range(symmetry_type):
        angle = step * angle_step
        target = _normalize(
            np.cos(angle) * seed_dir
            + np.sin(angle) * np.cross(axis, seed_dir)
            + (1.0 - np.cos(angle)) * np.dot(axis, seed_dir) * axis
        )
        scores = monomer_dirs @ target
        scores[~available] = -np.inf
        frame = int(np.argmax(scores))
        frames.append(frame)
        available[frame] = False

    universe.trajectory[0]
    return frames


def _pick_reference_atoms(universe, frames, roi_selection=None):
    def _pick_atom_from_monomer(frame):
        universe.trajectory[frame]
        atoms = universe.select_atoms(roi_selection) if roi_selection else None
        if atoms is None or len(atoms) == 0:
            atoms = universe.select_atoms("protein")

        ca = atoms.select_atoms("name CA")
        candidates = ca if len(ca) > 0 else atoms
        center = atoms.positions.mean(axis=0)
        distances = np.linalg.norm(candidates.positions - center, axis=1)
        return candidates[int(np.argmin(distances))].index

    ref_indices = [_pick_atom_from_monomer(frame) for frame in frames]
    universe.trajectory[0]
    return ref_indices


def select_2fold_pairs(universe, ref_idx, ref_frame, partner_frame,
                        ref_chain_id=None, partner_chain_id=None,
                        min_pair_span=5.0, max_pair_span=25.0):
    """
    Select two atom pairs (A, A') and (B, B') for 2-fold plane fitting.

    The reference atom A is identified by ``ref_idx`` in ``ref_frame``.
    Its partner A' is the same (resname, resid, name) atom in ``partner_frame``.
    A second atom B is chosen from the same local region as A, preferring
    nearby CA/heavy atoms within a useful distance span from A.

    Candidate selection priority:
      1. Same-residue heavy atoms (CA, N, C, O, CB, ...).
      2. Nearby CA atoms from neighbouring residues.
      3. Any atom whose distance from A is within [min_pair_span, max_pair_span].
    Candidates are rejected if no matching partner B' exists in ``partner_frame``.

    Parameters
    ----------
    universe : MDAnalysis.Universe
    ref_idx : int
        0-based global atom index of the primary reference atom A.
    ref_frame : int
        Trajectory frame (0-based) containing A.
    partner_frame : int
        Trajectory frame containing A'.
    ref_chain_id : str or None
        Optional original chain ID containing A and B in ``ref_frame``.
    partner_chain_id : str or None
        Optional original chain ID containing A' and B' in ``partner_frame``.
    min_pair_span : float
        Minimum allowed distance between A and B (Angstrom). Default 5.0.
    max_pair_span : float
        Maximum allowed distance between A and B (Angstrom). Default 25.0.

    Returns
    -------
    points : np.ndarray of shape (4, 3) or None
        [pointA, pointA_prime, pointB, pointB_prime] if a suitable B is found,
        otherwise None.
    diagnostics : dict
        Keys: 'b_atom_index', 'b_resname', 'b_resid', 'b_name',
              'pair_distance', 'candidate_count', 'error' (on failure).
    """
    # Get pointA and its atom identity
    universe.trajectory[ref_frame]
    ref_atom = universe.select_atoms(f"protein and index {ref_idx}")
    if len(ref_atom) == 0:
        return None, {"error": f"ref_idx {ref_idx} not found in frame {ref_frame}"}
    pointA = ref_atom.positions[0].copy()
    ref_resname = ref_atom.resnames[0]
    ref_resid = ref_atom.resids[0]
    ref_name = ref_atom.names[0]

    # Get pointA' in partner_frame
    universe.trajectory[partner_frame]
    partner_sel = (
        f"protein and resname {ref_resname} and resid {ref_resid} "
        f"and name {ref_name}"
    )
    if partner_chain_id is not None:
        partner_sel += f" and chainid {partner_chain_id}"
    partner_atoms = universe.select_atoms(partner_sel)
    if len(partner_atoms) == 0:
        universe.trajectory[0]
        return None, {"error": f"A' not found in frame {partner_frame}"}
    pointA_prime = partner_atoms.positions[0].copy()

    # Build candidate list for B in ref_frame.
    # Constrain to atoms within ±5 residues of ref_resid — this ensures B
    # belongs to the same local 2-fold patch as A, not a distant region.
    universe.trajectory[ref_frame]
    min_resid = ref_resid - 5
    max_resid = ref_resid + 5
    candidate_sel = f"protein and resid {min_resid}:{max_resid}"
    if ref_chain_id is not None:
        candidate_sel += f" and chainid {ref_chain_id}"
    candidate_pool = universe.select_atoms(candidate_sel)
    # If the residue window matches nothing (unusual), widen to all protein
    if len(candidate_pool) == 0:
        fallback_sel = "protein"
        if ref_chain_id is not None:
            fallback_sel += f" and chainid {ref_chain_id}"
        candidate_pool = universe.select_atoms(fallback_sel)

    candidates = []
    mid_span = (min_pair_span + max_pair_span) / 2.0

    for atom in candidate_pool:
        # Skip the reference atom itself
        if (atom.resname == ref_resname and atom.resid == ref_resid
                and atom.name == ref_name):
            continue

        b_pos = atom.position
        dist = np.linalg.norm(b_pos - pointA)
        if dist < min_pair_span or dist > max_pair_span:
            continue

        # Check B' exists in partner_frame
        universe.trajectory[partner_frame]
        b_prime_sel = (
            f"protein and resname {atom.resname} and resid {atom.resid} "
            f"and name {atom.name}"
        )
        if partner_chain_id is not None:
            b_prime_sel += f" and chainid {partner_chain_id}"
        b_prime_atoms = universe.select_atoms(b_prime_sel)
        if len(b_prime_atoms) == 0:
            universe.trajectory[ref_frame]
            continue

        # Score: prefer same-residue, CA > backbone heavy > other, prefer mid-span
        score = 0.0
        if atom.resname == ref_resname and atom.resid == ref_resid:
            score += 4.0
        if atom.name == "CA":
            score += 2.0
        elif atom.name in ("N", "C", "O", "CB"):
            score += 1.0
        score -= abs(dist - mid_span) / mid_span  # penalty for edge distances

        candidates.append({
            "score": score,
            "position": b_pos.copy(),
            "b_prime_pos": b_prime_atoms.positions[0].copy(),
            "index": atom.index,
            "resname": atom.resname,
            "resid": atom.resid,
            "name": atom.name,
            "distance": dist,
        })
        universe.trajectory[ref_frame]

    universe.trajectory[0]

    if len(candidates) == 0:
        return None, {
            "candidate_count": 0,
            "error": (f"No suitable B atom found with distance in "
                      f"[{min_pair_span:.0f}, {max_pair_span:.0f}]A "
                      f"and valid partner copy in frame {partner_frame}"),
        }

    candidates.sort(key=lambda c: c["score"], reverse=True)
    best = candidates[0]

    points = np.array([
        pointA,
        pointA_prime,
        best["position"],
        best["b_prime_pos"],
    ])

    diagnostics = {
        "b_atom_index": int(best["index"]),
        "b_resname": best["resname"],
        "b_resid": best["resid"],
        "b_name": best["name"],
        "pair_distance": best["distance"],
        "candidate_count": len(candidates),
        "ref_frame": int(ref_frame),
        "partner_frame": int(partner_frame),
        "ref_chain_id": ref_chain_id,
        "partner_chain_id": partner_chain_id,
    }

    return points, diagnostics


def compute_local_plane_diagnostics(points, detected_axis=None):
    """
    Compute quality metrics for a locally-fitted plane relative to a
    detected global symmetry axis.

    Parameters
    ----------
    points : np.ndarray of shape (N, 3)
        Points used in the local plane fit.
    detected_axis : np.ndarray of shape (3,) or None
        Reference axis direction (unit vector) from global symmetry detection.

    Returns
    -------
    diag : dict
        Keys: 'plane_normal', 'plane_rmsd', 'n_points',
              'alignment_angle_deg' (only if detected_axis is not None).
    """
    normal, rmsd, sv = fit_plane_normal(points, reference_direction=detected_axis)
    # non_collinearity ratio: close to 1 = well-spread, near 0 = nearly collinear
    non_collinearity = sv[1] / sv[0] if sv[0] > 1e-12 else 0.0
    diag = {
        "plane_normal": normal,
        "plane_rmsd": rmsd,
        "n_points": len(points),
        "singular_values": sv,
        "non_collinearity": non_collinearity,
    }
    if detected_axis is not None:
        dot = np.clip(np.abs(np.dot(normal, detected_axis)), -1.0, 1.0)
        diag["alignment_angle_deg"] = np.rad2deg(np.arccos(dot))
    return diag


def auto_detect_reference(universe, symmetry_type, axis_index=0,
                          roi_selection=None, roi_frame=0):
    """
    Automatically detect reference atom indices for a given symmetry type.

    This is the main entry point.  It:
      1. Computes monomer COMs (one per MODEL) and the capsid-wide centre.
      2. Direct spherical search finds the symmetry axis.
      3. Finds three monomers at the correct angular positions around the axis.
      4. Returns reference atom indices from those three monomers.

    When used with the existing ``sliceCapsid.py`` pipeline, the three
    returned indices can be used directly as *pointA*, *pointB*, *pointC*
    without the nearest-neighbour distance heuristic.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        The capsid structure (already loaded from a multi-model PDB file).
    symmetry_type : int
        Symmetry axis type: 2, 3, or 5.
    axis_index : int
        Which axis among the sorted candidates to use (0 = best).  If
        ``roi_selection`` is provided, this refers to the ROI-ranked list.
        Otherwise it refers to the global symmetry-ranked list.
    roi_selection : str or None
        Optional MDAnalysis selection describing the ROI on ``roi_frame``.
    roi_frame : int
        Frame containing the physical ROI copy described by ``roi_selection``.

    Returns
    -------
    ref_indices : list of int (length 3)
        0-based atom indices for *pointA*, *pointB*, *pointC*.
    axis_direction : np.ndarray of shape (3,)
        The detected symmetry axis direction (unit vector), for diagnostics.
    ref_frames : list of int (length 3)
        Trajectory frames for each reference atom.
    """
    capsid_center, monomer_coms, monomer_dirs, all_axes = _compute_axis_context(
        universe
    )

    # Step 3: Select the best axis of the requested type
    key = _axis_key(symmetry_type)
    candidates = all_axes[key]

    if len(candidates) == 0:
        raise RuntimeError(f"No {key} axes found — this should not happen.")

    if roi_selection:
        roi_center = compute_roi_center(
            universe, roi_selection, roi_frame=roi_frame
        )
        ranked_axes = _rank_axes_for_roi_from_context(
            candidates,
            monomer_coms,
            monomer_dirs,
            capsid_center,
            roi_center,
            symmetry_type,
        )
        if axis_index >= len(ranked_axes):
            raise ValueError(
                f"axis_index={axis_index} but only {len(ranked_axes)} {key} "
                f"axes are available. Valid range: 0–{len(ranked_axes) - 1}."
            )
        best_axis = ranked_axes[axis_index]["axis"]
    else:
        if axis_index >= len(candidates):
            raise ValueError(
                f"axis_index={axis_index} but only {len(candidates)} {key} axes available. "
                f"Valid range: 0–{len(candidates) - 1}."
            )
        best_axis = candidates[axis_index]

    # Step 4: Find three monomers at the correct angular positions around axis
    ref_frames = _find_reference_frames(
        monomer_dirs,
        best_axis,
        symmetry_type,
        prefer_positive_axis=bool(roi_selection),
    )

    # Step 5: Pick a representative atom from each monomer.
    ref_indices = _pick_reference_atoms(
        universe, ref_frames, roi_selection=roi_selection
    )

    # Return to frame 0
    universe.trajectory[0]

    return ref_indices, best_axis, ref_frames


# ---------------------------------------------------------------------------
# Convenience: list all available axes with scores (for diagnostics / user
# selection)
# ---------------------------------------------------------------------------

def list_axes(universe, symmetry_type, roi_selection=None, roi_frame=0):
    """
    Return all candidate axes of *symmetry_type* with quality scores.

    Parameters
    ----------
    universe : MDAnalysis.Universe
    symmetry_type : int
        2, 3, or 5.

    Returns
    -------
    results : list of dict
        Each dict has keys: 'axis' (np.ndarray), 'score' (float),
        'ref_index' (int), 'ref_frame' (int).
        Sorted by descending score.
    """
    capsid_center, monomer_coms, monomer_dirs, all_axes = _compute_axis_context(
        universe
    )
    key = _axis_key(symmetry_type)
    candidates = all_axes[key]

    results = []
    roi_rank_by_index = {}
    if roi_selection:
        roi_center = compute_roi_center(
            universe, roi_selection, roi_frame=roi_frame
        )
        roi_ranked = _rank_axes_for_roi_from_context(
            candidates,
            monomer_coms,
            monomer_dirs,
            capsid_center,
            roi_center,
            symmetry_type,
        )
        roi_rank_by_index = {
            item["axis_index"]: (rank, item)
            for rank, item in enumerate(roi_ranked)
        }
        if len(roi_rank_by_index) != len(candidates):
            raise RuntimeError("ROI ranking did not cover all candidate axes.")

    for i, ax in enumerate(candidates):
        if i in roi_rank_by_index:
            roi_rank, roi_item = roi_rank_by_index[i]
            ax = roi_item["axis"]
        else:
            roi_rank = None
            roi_item = {}

        score = _rotational_symmetry_score(ax, monomer_dirs, symmetry_type)

        ref_idx, ref_frame = find_reference_atom(
            universe, ax, capsid_center, symmetry_type
        )
        # Restore frame 0
        universe.trajectory[0]
        ref_atom = universe.select_atoms(f"protein and index {ref_idx}")
        chain_id = ref_atom.chainIDs[0] if len(ref_atom) > 0 else "?"

        results.append({
            "axis": ax,
            "score": score,
            "ref_index": ref_idx,
            "ref_frame": ref_frame,
            "chain_id": chain_id,
            "axis_index": i,
            "roi_rank": roi_rank,
            "roi_score": roi_item.get("roi_score"),
            "roi_angle_deg": roi_item.get("roi_angle_deg"),
            "roi_line_distance": roi_item.get("roi_line_distance"),
        })

    if roi_selection:
        results.sort(key=lambda x: x["roi_rank"])  # lower = closer ROI match
    else:
        results.sort(key=lambda x: x["score"])  # lower = better
    return results
