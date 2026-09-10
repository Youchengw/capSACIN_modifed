"""Exercise native startup and all three folds using the distributed engine."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


def main():
    app = Path(sys.argv[1]).resolve()
    executable_dir = app / "Contents/MacOS"
    subprocess.run(
        [str(executable_dir / "capsacin-studio")], check=True, timeout=30,
        env={**os.environ, "CAPSACIN_STARTUP_CHECK": "1"},
    )
    with tempfile.TemporaryDirectory(prefix="capsacin-release-smoke-") as workspace:
        request = {
            "id": "release-smoke", "operation": "prepare_all_previews",
            "params": {"input_path": str(app / "Contents/Resources/pdb/1k3v.pdb")},
        }
        start = time.monotonic()
        result = subprocess.run(
            [str(executable_dir / "capsacin-sidecar")],
            input=json.dumps(request) + "\n", text=True, capture_output=True,
            timeout=300, env={**os.environ, "CAPSACIN_WORKSPACE": workspace},
        )
        if result.returncode:
            raise RuntimeError(f"Packaged engine failed ({result.returncode}):\n{result.stderr}")
        messages = [json.loads(line) for line in result.stdout.splitlines() if line.strip()]
        ready = [item for item in messages if item.get("type") == "ready"]
        responses = [item for item in messages if item.get("type") == "result"]
        assert len(ready) == 1, f"Expected one server, got {len(ready)}"
        assert len(responses) == 1, f"Expected one result, got {len(responses)}"
        response = responses[0]
        assert response["id"] == request["id"] and response["status"] == "ok", response
        data = response["data"]
        assert not data["errors"], data["errors"]
        assert set(data["previews"]) == {"2", "3", "5"}, data["previews"].keys()
        for preview in data["previews"].values():
            cif = Path(preview["aligned_cif_path"])
            assert cif.is_file() and cif.stat().st_size > 0, cif
        print(f"Packaged startup and 2/3/5-fold previews passed ({time.monotonic() - start:.1f}s including engine startup)")


if __name__ == "__main__":
    main()
