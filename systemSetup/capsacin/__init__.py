"""capSACIN: Capsid Surface Abstraction and Computationally-Induced Nanofragmentation."""

from . import definePlane
from .definePlane import (
    fit_plane_normal,
    select_5fold_ring,
    compute_pentagon_diagnostics,
)

# Pipeline (Phase 1 refactoring)
from . import pipeline
from .pipeline import CapsidPipeline

# Protocol dataclasses
from . import protocol
from .protocol import (
    SliceRequest,
    LoadResult,
    AxisCandidate,
    AxisResult,
    PlaneResult,
    PDBInfoResult,
    AlignResult,
    FormattedResult,
    SliceDataResult,
    CleanupResult,
    SliceOutput,
    PipelineCancelledError,
)
