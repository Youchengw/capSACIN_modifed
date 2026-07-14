"""Protocol dataclasses for the capSACIN pipeline.

These define the structured data that flows between CapsidPipeline stages.
Internal dataclasses (containing MDAnalysis Universes or DataFrames) are
used within the pipeline; serializable dataclasses are used by the sidecar
and desktop app.
"""

from __future__ import annotations

import threading
from dataclasses import dataclass, field
from typing import Any


# ---------------------------------------------------------------------------
# Exception
# ---------------------------------------------------------------------------


class PipelineCancelledError(Exception):
    """Raised when a pipeline stage is cancelled via the cancel event."""


# ---------------------------------------------------------------------------
# Request (serializable)
# ---------------------------------------------------------------------------


@dataclass
class SliceRequest:
    """All parameters needed to run the full slice pipeline."""

    input_path: str
    workspace_path: str = ""
    symmetry: int = 5  # 2, 3, or 5
    weight: float = 0.5  # 0..1 slicing height fraction
    auto: bool = True
    axis_index: int = 0
    roi_selection: str | None = None
    roi_frame: int = 0
    ref_indices: list[int] | None = None
    legacy_plane_heuristic: bool = False
    plot: bool = False


# ---------------------------------------------------------------------------
# Internal pipeline dataclasses
# ---------------------------------------------------------------------------


@dataclass
class LoadResult:
    """Result of loading a PDB structure."""

    universe: Any  # MDAnalysis Universe (not serializable)
    orig_chains: list[str]
    n_monomers: int
    total_atoms: int
    total_residues: int


@dataclass
class AxisCandidate:
    """A single symmetry-axis candidate (serializable)."""

    axis: list[float]  # 3-element unit vector
    score: float  # symmetry score (lower = better)
    axis_index: int
    chain_id: str | None = None
    ref_frame: int | None = None
    ref_index: int | None = None
    # ROI-related fields (only populated when roi_selection is set)
    roi_score: float | None = None
    roi_angle_deg: float | None = None
    roi_line_distance: float | None = None


@dataclass
class AxisResult:
    """Result of symmetry-axis detection."""

    axis_dir: list[float]  # 3-element unit vector for alignment
    ref_indices: list[int]  # [pointA, pointB, pointC] 0-based atom indices
    ref_frames: list[int]  # trajectory frames for each reference atom
    candidates: list[AxisCandidate] = field(default_factory=list)
    selected_index: int = 0


@dataclass
class PlaneResult:
    """Result of plane-point selection and normal computation."""

    pointA: Any  # np.ndarray (3,) — primary reference point
    pointB: Any  # np.ndarray (3,) — secondary reference point
    pointC: Any  # np.ndarray (3,) — tertiary reference point
    plane_points: Any  # np.ndarray (N, 3) — points for SVD plane fitting
    normal_vector: Any  # np.ndarray (3,) — unit normal for alignment
    # Diagnostics
    plane_rmsd: float | None = None
    non_collinearity: float | None = None
    local_vs_global_angle_deg: float | None = None


@dataclass
class PDBInfoResult:
    """Extracted PDB atom/chain/residue information."""

    chains: Any  # pd.DataFrame with PDB columns plus viewer_chain_index
    coords: Any  # np.ndarray (N, 3)
    capsid_com: Any  # np.ndarray (3,)
    amino_acid_checks: dict[str, int]  # {resname: expected_atom_count}


@dataclass
class AlignResult:
    """Result of coordinate alignment."""

    aligned_coords: Any  # np.ndarray (N, 3) — rotated coordinates
    rotation_matrix: Any  # np.ndarray (3, 3)
    select_points: Any  # np.ndarray (M, 3) — aligned reference points


@dataclass
class FormattedResult:
    """Aligned, cleaned, and PDB-formatted chains."""

    chains: Any  # pd.DataFrame — formatted PDB-ready DataFrame
    aligned_coords_df: Any  # pd.DataFrame with x, y, z columns


@dataclass
class SliceDataResult:
    """Result of weight-based slicing."""

    sliced_chains: Any  # pd.DataFrame — chains above the z-cutoff
    cutoff: float
    original_lengths: list[int]  # expected residue count per chain


@dataclass
class CleanupResult:
    """Result of broken residue/chain cleanup."""

    cleaned_chains: Any  # pd.DataFrame
    fragmented_chains_dropped: int
    fragmented_residues_dropped: int
    original_lengths: list[int]


# ---------------------------------------------------------------------------
# Final output (serializable)
# ---------------------------------------------------------------------------


@dataclass
class SliceOutput:
    """Final serializable result of a slice operation."""

    output_pdb_path: str
    preview_cif_path: str = ""
    input_atoms: int = 0
    input_residues: int = 0
    input_chains: int = 0
    output_atoms: int = 0
    output_residues: int = 0
    output_chains: int = 0
    retained_percentage: float = 0.0
    axis_dir: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    diagnostics: dict = field(default_factory=dict)
    log_lines: list[str] = field(default_factory=list)
    # Internal viewer payload. It preserves pre-rename MODEL/chain identity and
    # is consumed by the desktop sidecar, not included in the JSON response.
    viewer_chains: Any = field(default=None, repr=False, compare=False)


# ---------------------------------------------------------------------------
# Sidecar protocol types (serializable)
# ---------------------------------------------------------------------------


@dataclass
class InspectResult:
    """Result of inspect_structure operation."""

    chains: list[dict] = field(default_factory=list)
    n_models: int = 0
    total_atoms: int = 0
    total_residues: int = 0
    valid_residue_ranges: dict[str, list[int]] = field(default_factory=dict)


@dataclass
class PreparePreviewResult:
    """Result of prepare_preview operation."""

    axis_candidates: list[AxisCandidate] = field(default_factory=list)
    selected_axis: AxisCandidate | None = None
    aligned_cif_path: str = ""
    diagnostics: dict = field(default_factory=dict)


@dataclass
class SidecarProgress:
    """A progress event emitted during sidecar operations."""

    id: str
    stage: str
    fraction: float
    message: str = ""
