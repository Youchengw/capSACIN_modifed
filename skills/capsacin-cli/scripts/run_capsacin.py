#!/usr/bin/env python3
"""Launch the existing capSACIN CLI with portable input and working paths."""

import argparse
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile


def find_repo(explicit):
    if explicit:
        candidates = [Path(explicit).expanduser().resolve()]
    else:
        cwd = Path.cwd()
        candidates = [cwd, *cwd.parents, *Path(__file__).resolve().parents]
    for candidate in candidates:
        if (candidate / "systemSetup" / "sliceCapsid.py").is_file():
            return candidate
    raise ValueError("Cannot locate capSACIN; set --repo or CAPSACIN_REPO to its checkout.")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog="Put wrapper options before the command, then -- and native CLI options.",
    )
    parser.add_argument("--repo", default=os.environ.get("CAPSACIN_REPO"))
    parser.add_argument(
        "--python", default=os.environ.get("CAPSACIN_PYTHON", sys.executable),
        help="Python executable with capSACIN dependencies (default: this interpreter)",
    )
    parser.add_argument("--input", required=True, help="Path to a local capsid PDB")
    parser.add_argument(
        "--work-dir", help="New directory for a slice/axis run; default: unique directory in cwd",
    )
    parser.add_argument("command", choices=("inspect", "slice"))
    parser.add_argument("cli_args", nargs=argparse.REMAINDER, help="Native CLI arguments")
    args = parser.parse_args(argv)
    native_args = args.cli_args[1:] if args.cli_args[:1] == ["--"] else args.cli_args

    try:
        repo = find_repo(args.repo)
        source = Path(args.input).expanduser().resolve()
        if not source.is_file():
            raise ValueError(f"Input file not found: {source}")
        python = shutil.which(os.path.expanduser(args.python))
        if python is None:
            raise ValueError(f"Python executable not found: {args.python}")
        # Keep a virtual environment's executable path; resolving its symlink
        # can select the base interpreter and lose that environment's packages.
        python = os.path.abspath(python)
        setup = repo / "systemSetup"
        if args.command == "inspect":
            if args.work_dir:
                raise ValueError("inspect does not use --work-dir; remove that option.")
            command = [python, "-u", str(setup / "inspectResidueRanges.py"),
                       str(source), *native_args]
            return subprocess.run(command, check=False).returncode

        if args.work_dir:
            work_dir = Path(args.work_dir).expanduser().resolve()
            work_dir.mkdir(parents=True, exist_ok=False)
        else:
            work_dir = Path(tempfile.mkdtemp(prefix="capsacin-", dir=Path.cwd()))
        (work_dir / "input").mkdir()
        (work_dir / "output").mkdir()
        staged = work_dir / "input" / f"{source.stem}.pdb"
        shutil.copyfile(source, staged)

        # --pdb is managed here because the native CLI requires input/<prefix>.pdb.
        command = [python, "-u", str(setup / "sliceCapsid.py"), *native_args,
                   f"--pdb={source.stem}"]
        env = os.environ.copy()
        env.setdefault("MPLBACKEND", "Agg")
        env.setdefault("MPLCONFIGDIR", str(work_dir / ".matplotlib"))
        log_path = work_dir / "run.log"
        with log_path.open("w", encoding="utf-8", buffering=1) as log:
            def emit(message):
                print(message, flush=True)
                log.write(message + "\n")

            emit(f"Input: {source}")
            emit(f"Repository: {repo}")
            emit(f"Working directory: {work_dir}")
            emit(f"Command: {shlex.join(command)}")
            with subprocess.Popen(
                command, cwd=work_dir, env=env,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding="utf-8", errors="replace", bufsize=1,
            ) as process:
                for line in process.stdout:
                    print(line, end="", flush=True)
                    log.write(line)
                returncode = process.wait()
            emit(f"Exit status: {returncode}")
            emit(f"Log: {log_path}")
            for artifact in sorted((work_dir / "output").glob("*.pdb")):
                emit(f"Output PDB: {artifact} ({artifact.stat().st_size} bytes)")
        return returncode
    except (OSError, ValueError) as exc:
        print(f"capSACIN launcher: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
