"""Batch preparation preserves independent alignment and isolates fold failures."""
import threading
from pathlib import Path

import numpy as np
import pytest

from capsacin import findSymmetryAxes as axes
from capsacin.pipeline import CapsidPipeline
from capsacin.protocol import PipelineCancelledError
from sidecar.handlers import handle_prepare_all_previews


@pytest.mark.parametrize("roi", [None, "protein"])
def test_batch_matches_independent_alignment(pdb_9jjh, tmp_path, monkeypatch, roi):
    calls = {"load_structure": 0, "extract_pdb_info": 0, "_compute_axis_context": 0}
    for owner, name in [(CapsidPipeline, "load_structure"),
                        (CapsidPipeline, "extract_pdb_info"), (axes, "_compute_axis_context")]:
        original = getattr(owner, name)

        def counted(*args, _name=name, _original=original, **kwargs):
            calls[_name] += 1
            return _original(*args, **kwargs)

        monkeypatch.setattr(owner, name, counted)
    # Preview menus must not trigger the costly CLI reference diagnostics.
    def unexpected_reference(*args, **kwargs):
        raise AssertionError("Preview calculated an unused candidate reference")
    monkeypatch.setattr(axes, "find_reference_atom", unexpected_reference)
    progress = []
    result = handle_prepare_all_previews(
        {"input_path": pdb_9jjh, "symmetry": 5, "roi_selection": roi,
         "axis_indices": {"2": 1, "3": 0, "5": 2}},
        tmp_path, lambda stage, fraction, message="": progress.append(fraction),
    )
    assert result["errors"] == {}
    assert set(result["previews"]) == {"2", "3", "5"}
    assert calls == {"load_structure": 1, "extract_pdb_info": 1, "_compute_axis_context": 1}
    assert progress == sorted(progress)
    assert progress[-1] == 1
    paths = [item["aligned_cif_path"] for item in result["previews"].values()]
    assert len(set(paths)) == 3

    for fold, index in [(2, 1), (3, 0), (5, 2)]:
        pipeline = CapsidPipeline()
        pipeline.load_structure(pdb_9jjh)
        axis = pipeline.detect_axis(fold, axis_index=index, roi_selection=roi)
        plane = pipeline.select_plane_points(axis, fold, auto_mode=True)
        pdb_info = pipeline.extract_pdb_info()
        aligned = pipeline.align_coordinates(pdb_info, plane, fold, auto_mode=True)
        formatted = pipeline.format_aligned_pdb(pdb_info, aligned)
        actual = result["previews"][str(fold)]
        np.testing.assert_array_equal(actual["selected_axis"], np.round(axis.axis_dir, 6))
        coords = np.asarray([
            line.split()[12:15]
            for line in Path(actual["aligned_cif_path"]).read_text().splitlines()
            if line.startswith("ATOM ")
        ], dtype=float)
        np.testing.assert_array_equal(coords, formatted.chains[["x", "y", "z"]].astype(float))
        expected_frames = axes.find_symmetry_related_frames(
            pipeline._u, plane.normal_vector, fold, seed_frame=0 if roi else None,
        )
        assert actual["diagnostics"]["roi_symmetry_frames"] == expected_frames


def test_one_invalid_candidate_keeps_other_folds(pdb_9jjh, tmp_path):
    result = handle_prepare_all_previews(
        {"input_path": pdb_9jjh, "axis_indices": {"3": 999}}, tmp_path,
    )
    assert set(result["previews"]) == {"2", "5"}
    assert set(result["errors"]) == {"3"}
    assert "axis_index=999" in result["errors"]["3"]


def test_cancel_stops_before_preparing_the_next_fold(pdb_9jjh, tmp_path):
    cancel = threading.Event()

    def progress(stage, fraction, message=""):
        if message == "Preview ready":
            cancel.set()

    with pytest.raises(PipelineCancelledError):
        handle_prepare_all_previews({"input_path": pdb_9jjh}, tmp_path, progress, cancel)
    assert len(list(tmp_path.glob("prepare_*/aligned_preview.cif"))) == 1
