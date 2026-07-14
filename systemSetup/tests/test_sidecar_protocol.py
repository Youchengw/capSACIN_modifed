"""Focused regression tests for the desktop sidecar protocol."""

from __future__ import annotations

import io
import json
import sys
import time

import pandas as pd

from sidecar.main import SidecarServer
from sidecar.handlers import _get_workspace
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
