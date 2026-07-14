"""Unit and integration tests for capsacin.pipeline.CapsidPipeline."""
import os
import threading

import numpy as np
import pytest

from capsacin.pipeline import CapsidPipeline
from capsacin import findSymmetryAxes
from capsacin.protocol import (
    PipelineCancelledError,
    SliceRequest,
)


class TestCapsidPipeline:
    """Tests for the CapsidPipeline class."""

    def test_create_pipeline(self):
        """Pipeline can be instantiated without arguments."""
        p = CapsidPipeline()
        assert p is not None

    def test_load_small_pdb(self, pdb_9jjh):
        """Loading a valid PDB returns a LoadResult with expected counts."""
        p = CapsidPipeline()
        result = p.load_structure(pdb_9jjh)

        assert result.n_monomers > 0
        assert result.total_atoms > 0
        assert result.total_residues > 0
        assert len(result.orig_chains) == result.n_monomers
        # 9jjh has 60 frames, 1 chain per frame (NM_001), 1298 atoms per chain
        assert result.total_atoms == 77880
        assert result.n_monomers == 1

    def test_load_large_pdb(self, pdb_1k3v):
        """Loading 1k3v returns valid chain counts."""
        p = CapsidPipeline()
        result = p.load_structure(pdb_1k3v)

        assert result.n_monomers > 0
        assert result.total_atoms > 100000
        # 1k3v has 258,360 atoms (60 frames x 2 chains x ~2153 atoms)
        assert result.total_atoms == 258360

    def test_detect_axis_auto_5fold(self, pdb_9jjh):
        """Auto-detection finds 5-fold axes for 9jjh."""
        p = CapsidPipeline()
        p.load_structure(pdb_9jjh)

        result = p.detect_axis(symmetry=5, auto=True)

        assert result.axis_dir is not None
        assert len(result.axis_dir) == 3
        # Check it's approximately a unit vector
        norm = np.linalg.norm(result.axis_dir)
        assert abs(norm - 1.0) < 0.01
        assert len(result.ref_indices) == 3
        assert len(result.ref_frames) == 3

    def test_complete_5fold_roi_frames_include_seed(self, pdb_9jjh):
        """Viewer ROI expansion finds five unique symmetry-related MODELs."""
        p = CapsidPipeline()
        p.load_structure(pdb_9jjh)
        axis = p.detect_axis(symmetry=5, auto=True).axis_dir

        frames = findSymmetryAxes.find_symmetry_related_frames(
            p._u, axis, symmetry_type=5, seed_frame=0
        )

        assert len(frames) == 5
        assert len(set(frames)) == 5
        assert frames[0] == 0

    def test_global_3fold_roi_frames_follow_selected_axis(self, pdb_1k3v):
        """Global ROI display derives its three copies from the selected axis."""
        p = CapsidPipeline()
        p.load_structure(pdb_1k3v)
        axis = np.asarray(p.detect_axis(symmetry=3, auto=True).axis_dir)

        frames = findSymmetryAxes.find_symmetry_related_frames(
            p._u, axis, symmetry_type=3
        )
        center = findSymmetryAxes.compute_capsid_center(p._u)
        coms = findSymmetryAxes.compute_monomer_coms(p._u)
        dirs = coms - center
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
        closest_frame = int(np.argmax(dirs @ axis))

        assert len(frames) == 3
        assert len(set(frames)) == 3
        assert frames[0] == closest_frame
        assert np.all(dirs[frames] @ axis > 0.9)

    def test_detect_axis_auto_2fold(self, pdb_1k3v):
        """Auto-detection finds 2-fold axes for 1k3v."""
        p = CapsidPipeline()
        p.load_structure(pdb_1k3v)

        result = p.detect_axis(symmetry=2, auto=True)

        assert len(result.axis_dir) == 3
        assert len(result.ref_indices) == 3

    def test_detect_axis_auto_3fold(self, pdb_1k3v):
        """Auto-detection finds 3-fold axes for 1k3v."""
        p = CapsidPipeline()
        p.load_structure(pdb_1k3v)

        result = p.detect_axis(symmetry=3, auto=True)

        assert len(result.axis_dir) == 3
        assert len(result.ref_indices) == 3

    def test_detect_axis_manual_requires_refindex(self, pdb_9jjh):
        """Manual mode raises if no ref_indices provided."""
        p = CapsidPipeline()
        p.load_structure(pdb_9jjh)

        with pytest.raises(ValueError, match="Manual mode requires"):
            p.detect_axis(symmetry=5, auto=False)

    def test_detect_axis_list_axes(self, pdb_9jjh):
        """list_axes returns candidate list."""
        p = CapsidPipeline()
        p.load_structure(pdb_9jjh)

        result = p.detect_axis(symmetry=5, auto=False, list_axes=True)

        assert len(result.candidates) > 0

    def test_invalid_symmetry_raises(self, pdb_9jjh):
        """Unsupported symmetry raises ValueError."""
        p = CapsidPipeline()
        p.load_structure(pdb_9jjh)

        with pytest.raises(ValueError, match="Unsupported symmetry"):
            p.detect_axis(symmetry=4, auto=False, list_axes=True)

    def test_extract_pdb_info(self, pdb_9jjh):
        """Extract PDB info produces correct DataFrame columns."""
        p = CapsidPipeline()
        p.load_structure(pdb_9jjh)

        result = p.extract_pdb_info()

        expected_cols = {"atom", "idx", "name", "resname", "chain",
                         "resids", "x", "y", "z", "occ", "b", "type",
                         "viewer_chain_index"}
        assert set(result.chains.columns) == expected_cols

        # Coordinate count should match atom count
        assert result.coords.shape[0] == p._load_result.total_atoms
        assert result.coords.shape[1] == 3

        # amino_acid_checks should have entries for common residues
        assert "ALA" in result.amino_acid_checks
        assert len(result.amino_acid_checks) > 0

    def test_full_pipeline_9jjh(self, pdb_9jjh, workspace):
        """Full pipeline on 9jjh produces valid output."""
        request = SliceRequest(
            input_path=pdb_9jjh,
            workspace_path=workspace,
            symmetry=5,
            weight=0.5,
            auto=True,
        )

        pipeline = CapsidPipeline()
        result = pipeline.run_full(request)

        assert result.output_atoms > 0
        assert result.output_atoms < result.input_atoms
        assert 0.0 < result.retained_percentage < 100.0
        assert result.output_chains > 0
        assert os.path.exists(result.output_pdb_path)

        # Verify PDB is readable by checking atom count in file
        with open(result.output_pdb_path) as f:
            atom_lines = [l for l in f if l.startswith("ATOM")]
        assert len(atom_lines) == result.output_atoms

    def test_full_pipeline_weight_extremes(self, pdb_9jjh, workspace):
        """Weight=0 keeps everything, weight=1 removes everything."""
        # weight = 0: keep everything
        req0 = SliceRequest(
            input_path=pdb_9jjh, workspace_path=workspace,
            symmetry=5, weight=0.0, auto=True,
        )
        r0 = CapsidPipeline().run_full(req0)
        assert r0.retained_percentage > 95.0  # nearly everything kept

        # weight = 1: remove everything
        req1 = SliceRequest(
            input_path=pdb_9jjh, workspace_path=workspace,
            symmetry=5, weight=1.0, auto=True,
        )
        r1 = CapsidPipeline().run_full(req1)
        assert r1.retained_percentage < 10.0  # nearly everything removed

    def test_full_pipeline_1k3v_all_folds(self, pdb_1k3v, workspace):
        """1k3v works with all three symmetry folds."""
        for sym in [2, 3, 5]:
            request = SliceRequest(
                input_path=pdb_1k3v, workspace_path=workspace,
                symmetry=sym, weight=0.55, auto=True,
            )
            result = CapsidPipeline().run_full(request)
            assert result.output_atoms > 0
            assert result.output_chains > 0

    def test_cancellation(self, pdb_9jjh):
        """Setting cancel_event stops the pipeline."""
        cancel = threading.Event()
        pipeline = CapsidPipeline(cancel_event=cancel)
        pipeline.load_structure(pdb_9jjh)

        # Cancel before a long operation
        cancel.set()

        with pytest.raises(PipelineCancelledError):
            pipeline.extract_pdb_info()

    def test_progress_callback(self, pdb_9jjh):
        """Progress callback receives stage updates."""
        stages_seen = []

        def track(stage, fraction):
            stages_seen.append((stage, fraction))

        pipeline = CapsidPipeline(progress_callback=track)
        pipeline.load_structure(pdb_9jjh)
        pipeline.detect_axis(symmetry=5, auto=True)

        assert len(stages_seen) > 0
        # All fractions should be in [0, 1]
        for _, frac in stages_seen:
            assert 0.0 <= frac <= 1.0

    def test_run_cli_preserves_output_naming(self, pdb_9jjh):
        """run_cli produces output with legacy naming convention."""
        import argparse

        args = argparse.Namespace(
            pdb="9jjh",
            symmetry=5,
            weight=0.5,
            auto=True,
            axis_index=0,
            roi_selection=None,
            roi_frame=0,
            refindex=[1058],
            legacy_plane_heuristic=False,
            plot=False,
            list_axes=False,
        )

        # run_cli is a static method
        exit_code = CapsidPipeline.run_cli(args)
        assert exit_code == 0

        expected_path = "output/9jjh-sliced-sym5-w0.5.pdb"
        assert os.path.exists(expected_path)
