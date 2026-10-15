"""Focused regression tests for the desktop sidecar protocol."""

from __future__ import annotations

import io
import json
import sys
import time
from pathlib import Path

import pandas as pd

from sidecar.main import SidecarServer
from sidecar.handlers import (
    _get_workspace, _chain_catalog, _write_chain_selection_source,
    handle_apply_chain_selection, handle_inspect_structure,
)
from sidecar.mmcif_writer import write_aligned_mmcif, write_sliced_mmcif


def test_concurrent_requests_keep_their_original_ids(monkeypatch, tmp_path):
    """A later stdin line must not overwrite an earlier worker's request id."""
    monkeypatch.setenv("CAPSACIN_WORKSPACE", str(tmp_path))
    server = SidecarServer()

    def fake_dispatch(request_id, operation, params, cancel_event):
        time.sleep(params["delay"])
        return {"request_id": request_id, "operation": operation}

    monkeypatch.setattr(server, "_dispatch", fake_dispatch)
    requests = [
        {"id": "slow", "operation": "first", "params": {"delay": 0.03}},
        {"id": "fast", "operation": "second", "params": {"delay": 0.01}},
        {"id": "stop", "operation": "shutdown", "params": {}},
    ]
    stdin = io.StringIO("".join(json.dumps(item) + "\n" for item in requests))
    stdout = io.StringIO()
    monkeypatch.setattr(sys, "stdin", stdin)
    monkeypatch.setattr(sys, "stdout", stdout)

    server.run()

    messages = [json.loads(line) for line in stdout.getvalue().splitlines()]
    results = {
        message["id"]: message
        for message in messages
        if message.get("type") == "result"
    }
    assert results["slow"]["data"] == {
        "request_id": "slow",
        "operation": "first",
    }
    assert results["fast"]["data"] == {
        "request_id": "fast",
        "operation": "second",
    }
    assert results["stop"]["status"] == "ok"


def test_unknown_operation_returns_correlated_error(monkeypatch, tmp_path):
    monkeypatch.setenv("CAPSACIN_WORKSPACE", str(tmp_path))
    server = SidecarServer()
    stdin = io.StringIO(
        json.dumps({"id": "bad-1", "operation": "does_not_exist", "params": {}})
        + "\n"
        + json.dumps({"id": "stop", "operation": "shutdown", "params": {}})
        + "\n"
    )
    stdout = io.StringIO()
    monkeypatch.setattr(sys, "stdin", stdin)
    monkeypatch.setattr(sys, "stdout", stdout)

    server.run()

    messages = [json.loads(line) for line in stdout.getvalue().splitlines()]
    error = next(message for message in messages if message.get("id") == "bad-1")
    assert error["status"] == "error"
    assert "Unknown operation" in error["error"]["message"]


def test_sliced_viewer_preserves_model_and_chain_identity(tmp_path):
    """ROI selection must still address M{frame}_{chain} after slicing."""
    rows = pd.DataFrame([
        {
            "atom": "ATOM", "idx": 1, "type": "C", "name": "CA",
            "resname": "ALA", "chain": "A", "resids": 42,
            "x": 1.0, "y": 2.0, "z": 3.0, "occ": 1.0,
        },
        {
            "atom": "ATOM", "idx": 2, "type": "C", "name": "CA",
            "resname": "ALA", "chain": "C", "resids": 42,
            "x": 4.0, "y": 5.0, "z": 6.0, "occ": 1.0,
        },
    ])
    output = tmp_path / "sliced.cif"

    write_sliced_mmcif(rows, str(output), orig_chains=["A", "B"])

    text = output.read_text()
    assert " M0_A 1 42 " in text
    assert " M1_A 1 42 " in text


def test_viewer_workspaces_are_unique_across_sequential_requests(tmp_path):
    """Mol* cache keys must never point at files overwritten by a later job."""
    first = _get_workspace(tmp_path, "prepare")
    second = _get_workspace(tmp_path, "prepare")

    assert first != second
    assert first.is_dir()
    assert second.is_dir()


def test_structure_inspection_counts_glup_as_protein(pdb_glup):
    """Desktop inspection must report GLUP atoms and residue 567."""
    result = handle_inspect_structure({"input_path": pdb_glup})

    assert result["total_atoms"] == 6
    assert result["total_residues"] == 2
    assert result["chains"] == [{
        "id": "A",
        "residue_count": 2,
        "min_resid": 1,
        "max_resid": 567,
        "atom_count": 6,
    }]


def test_aligned_viewer_uses_numeric_index_for_unicode_whitespace_chain(tmp_path):
    """8des-like flat chains must not disappear when chr(65+n) is whitespace."""
    rows = pd.DataFrame([{
        "atom": "ATOM", "idx": 1, "type": "C", "name": "CA",
        "resname": "ALA", "chain": "\x85", "resids": 42,
        "x": 1.0, "y": 2.0, "z": 3.0, "occ": 1.0, "b": 0.0,
        "viewer_chain_index": 68,
    }])
    output = tmp_path / "aligned.cif"

    write_aligned_mmcif(
        rows,
        orig_chains=["A", "E"],
        n_frames=60,
        output_path=str(output),
    )

    # flat index 68 = frame 34, first original chain A
    assert " M34_A 1 42 " in output.read_text()


def test_chain_catalog_and_manual_export_use_complete_aligned_chains(tmp_path):
    rows = pd.DataFrame([
        {"atom": "ATOM", "idx": "1", "type": "C", "name": "CA", "resname": "ALA",
         "chain": "A", "resids": "1", "x": "0.000", "y": "0.000", "z": "5.000",
         "occ": "1.00", "b": "0.00", "viewer_chain_index": "0"},
        {"atom": "ATOM", "idx": "2", "type": "C", "name": "CB", "resname": "ALA",
         "chain": "A", "resids": "1", "x": "1.000", "y": "0.000", "z": "7.000",
         "occ": "1.00", "b": "0.00", "viewer_chain_index": "0"},
        {"atom": "ATOM", "idx": "3", "type": "C", "name": "CA", "resname": "ALA",
         "chain": "B", "resids": "1", "x": "0.000", "y": "1.000", "z": "25.000",
         "occ": "1.00", "b": "0.00", "viewer_chain_index": "1"},
        {"atom": "ATOM", "idx": "4", "type": "C", "name": "CB", "resname": "ALA",
         "chain": "B", "resids": "1", "x": "1.000", "y": "1.000", "z": "27.000",
         "occ": "1.00", "b": "0.00", "viewer_chain_index": "1"},
    ])
    catalog = _chain_catalog(rows, ["X"], selected_indices=[1], weight_selected_indices=[1])
    assert [item["id"] for item in catalog] == ["M1_X", "M0_X"]
    assert [item["z_mean"] for item in catalog] == [26.0, 6.0]
    assert catalog[0]["selected"] is True
    assert catalog[1]["selected"] is False

    source = tmp_path / "chains.pkl"
    _write_chain_selection_source(source, rows, ["X"], {
        "input_atoms": 4, "input_residues": 2, "input_chains": 2,
        "axis": [0, 0, 1], "diagnostics": {"cutoff": 20}, "warnings": [],
        "weight_selected_indices": [1], "source_stem": "capsid",
    })
    result = handle_apply_chain_selection({
        "chain_source_path": str(source), "selected_chain_indices": [0],
    }, tmp_path)
    assert result["statistics"]["output_atoms"] == 2
    assert result["statistics"]["output_chains"] == 1
    assert result["chain_candidates"][0]["selected"] is False
    assert result["chain_candidates"][1]["selected"] is True
    assert " M0_X 1 1 " in Path(result["sliced_cif_path"]).read_text()
    assert " M1_X 1 1 " not in Path(result["sliced_cif_path"]).read_text()
