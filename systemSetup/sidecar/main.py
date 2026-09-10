#!/usr/bin/env python3
"""capSACIN sidecar — JSON-lines protocol process.

Communicates with the desktop app via stdin/stdout using one JSON object per line.

Protocol:
  Input  (stdin):  {"id":"...","operation":"...","params":{...}}
  Output (stdout): {"type":"progress","id":"...","stage":"...","fraction":0.5,"message":"..."}
  Output (stdout): {"type":"result","id":"...","status":"ok","data":{...}}
  Output (stdout): {"type":"result","id":"...","status":"error","error":{"code":"...","message":"..."}}

Operations:
  - inspect_structure   Read-only structural metadata
  - prepare_preview     Load, detect axis, align, generate aligned viewer mmCIF
  - run_slice           Full pipeline: prepare + slice + cleanup + write
  - cancel              Cancel a running operation
  - shutdown            Clean exit

Run directly:  python -m sidecar
"""
from __future__ import annotations

import json
import os
import sys
import threading
import traceback
from pathlib import Path
from typing import Any

from capsacin.pipeline import CapsidPipeline
from capsacin.protocol import PipelineCancelledError, SliceRequest

from .handlers import (
    handle_inspect_structure,
    handle_prepare_preview,
    handle_prepare_all_previews,
    handle_run_slice,
)


_STDOUT_LOCK = threading.Lock()


def _write_line(obj: dict) -> None:
    """Write a JSON object as a single line to stdout, flushed immediately."""
    with _STDOUT_LOCK:
        sys.stdout.write(json.dumps(obj, ensure_ascii=True, default=str) + "\n")
        sys.stdout.flush()


def _emit_progress(request_id: str, stage: str, fraction: float, message: str = "") -> None:
    _write_line({
        "type": "progress",
        "id": request_id,
        "stage": stage,
        "fraction": fraction,
        "message": message,
    })


def _emit_result(request_id: str, status: str, data: Any = None,
                 error: dict | None = None) -> None:
    payload: dict = {"type": "result", "id": request_id, "status": status}
    if data is not None:
        payload["data"] = data
    if error is not None:
        payload["error"] = error
    _write_line(payload)


def _emit_error(request_id: str, code: str, message: str) -> None:
    _emit_result(request_id, "error", error={"code": code, "message": message})


class SidecarServer:
    """Main sidecar event loop."""

    def __init__(self):
        self._cancel_events: dict[str, threading.Event] = {}
        self._workspace_base = Path(os.environ.get(
            "CAPSACIN_WORKSPACE",
            Path.home() / ".capsacin" / "workspaces",
        ))
        self._workspace_base.mkdir(parents=True, exist_ok=True)

    def run(self) -> None:
        """Block on stdin, dispatch JSON requests, write JSON responses to stdout."""
        _write_line({"type": "ready", "version": "1.0.0"})
        self._threads: list[threading.Thread] = []
        self._stdin_closed = False

        for line in sys.stdin:
            line = line.strip()
            if not line:
                continue

            try:
                request = json.loads(line)
            except json.JSONDecodeError as e:
                _emit_error("unknown", "INVALID_JSON", str(e))
                continue

            req_id = request.get("id", "unknown")
            operation = request.get("operation", "")
            params = request.get("params", {})

            if operation == "shutdown":
                _emit_result(req_id, "ok", data={"message": "Goodbye."})
                break

            elif operation == "cancel":
                target_id = params.get("target_id", req_id)
                if target_id in self._cancel_events:
                    self._cancel_events[target_id].set()
                    _emit_result(req_id, "ok", data={"cancelled": target_id})
                else:
                    _emit_error(req_id, "NOT_FOUND",
                                f"No active operation with id {target_id}")
                continue

            # Clean up completed threads
            self._threads = [t for t in self._threads if t.is_alive()]

            # Dispatch in a thread so we can handle cancellation
            cancel_event = threading.Event()
            self._cancel_events[req_id] = cancel_event

            def _runner(
                request_id=req_id,
                request_operation=operation,
                request_params=params,
                request_cancel_event=cancel_event,
            ):
                try:
                    result = self._dispatch(
                        request_id,
                        request_operation,
                        request_params,
                        request_cancel_event,
                    )
                    _emit_result(request_id, "ok", data=result)
                except PipelineCancelledError:
                    _emit_result(
                        request_id,
                        "cancelled",
                        data={"message": "Operation cancelled."},
                    )
                except Exception as e:
                    _emit_error(request_id, type(e).__name__.upper(), str(e))
                    traceback.print_exc(file=sys.stderr)
                finally:
                    self._cancel_events.pop(request_id, None)

            thread = threading.Thread(target=_runner, daemon=False)
            thread.start()
            self._threads.append(thread)

        # Wait for all pending operations to complete before exiting
        for thread in self._threads:
            thread.join(timeout=300)  # 5-minute max wait per thread

    def _dispatch(self, req_id: str, operation: str, params: dict,
                  cancel_event: threading.Event) -> dict:
        """Route to the appropriate handler."""
        if operation == "inspect_structure":
            return handle_inspect_structure(
                params,
                progress=lambda stage, frac, msg="": _emit_progress(req_id, stage, frac, msg),
                cancel=cancel_event,
            )
        elif operation in ("prepare_preview", "prepare_all_previews"):
            handler = (handle_prepare_all_previews if operation == "prepare_all_previews"
                       else handle_prepare_preview)
            return handler(
                params,
                workspace_base=self._workspace_base,
                progress=lambda stage, frac, msg="": _emit_progress(req_id, stage, frac, msg),
                cancel=cancel_event,
            )
        elif operation == "run_slice":
            return handle_run_slice(
                params,
                workspace_base=self._workspace_base,
                progress=lambda stage, frac, msg="": _emit_progress(req_id, stage, frac, msg),
                cancel=cancel_event,
            )
        else:
            raise ValueError(f"Unknown operation: {operation}")


def main():
    """Entry point for `python -m sidecar`."""
    server = SidecarServer()
    server.run()


if __name__ == "__main__":
    main()
