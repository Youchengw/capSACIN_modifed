---
name: capsacin-cli
description: "Run the capSACIN CLI to inspect capsid PDBs, list 2/3/5-fold axes, and slice capsids using automatic, ROI-aware, or manual reference selection. Use for capSACIN command-line execution and input troubleshooting, not general MD simulation or desktop UI control."
---

# capSACIN CLI

Use the existing capSACIN Python CLI to produce actual structures and report the
result. This skill contains instructions and a launcher; it requires a source
checkout and the capSACIN Python environment. It does not bundle the engine.

## Locate the engine

- Find the checkout containing `systemSetup/sliceCapsid.py`. The launcher accepts
  `--repo`, then `CAPSACIN_REPO`, or discovers the checkout above the current
  directory or this skill. Prefer an explicit path when multiple checkouts exist.
- Use the Python interpreter from the user's working capSACIN environment.
  `--python` or `CAPSACIN_PYTHON` selects it; otherwise the launcher uses its own
  interpreter. A copied skill does not imply that dependencies are installed.
- For first-time setup, read [references/cli.md](references/cli.md#environment).
  Reuse a working environment instead of rebuilding it on every invocation.

## Interpret the request

- Input: a local PDB containing spatial capsid copies as separate MODELs with
  compatible atom identities and ordering. These are symmetry copies, not MD
  snapshots. A single asymmetric unit, a single-MODEL full assembly, ordinary
  GRO, or mmCIF is not an equivalent input to this CLI. Renaming extensions does
  not convert structures. Bundled PDBs are useful for a requested demonstration.
- Use the requested symmetry (2, 3, or 5), slicing weight, and selection mode.
  Ask for missing scientific choices when context does not establish them;
  do not quietly use the native example defaults for a research result.
- Preserve explicit manual `--refindex` values. They are zero-based atom indices
  within a MODEL, not PDB serial numbers or residue IDs. Do not add `--auto` to
  a manual request. Use `--auto` for requested automatic or ROI-aware selection.
- ROI ranking requires a user-defined MDAnalysis selection and MODEL/frame.
  It finds an axis near that ROI; it does not discover a biological ROI.
  `--roi-frame` is the zero-based MODEL position in file order. The inspector's
  printed `MODEL - 1` mapping assumes consecutively numbered MODELs starting at 1.

## Execute

Use `scripts/run_capsacin.py` relative to this skill directory. Resolve its
absolute path before changing directories. Wrapper options go before `inspect`
or `slice`; arguments after `--` go directly to the native CLI.

```bash
# Replace paths with this installation and the requested input.
CAPSACIN_RUNNER=/path/to/capsacin-cli/scripts/run_capsacin.py
CAPSACIN_REPO=/path/to/capSACIN
CAPSACIN_PYTHON=/path/to/capSACIN-environment/bin/python
CAPSACIN_INPUT=/path/to/capsid.pdb

"$CAPSACIN_PYTHON" "$CAPSACIN_RUNNER" --repo "$CAPSACIN_REPO" \
  --input "$CAPSACIN_INPUT" inspect

# Example values only: automatic fivefold slice, weight 0.7.
"$CAPSACIN_PYTHON" "$CAPSACIN_RUNNER" --repo "$CAPSACIN_REPO" \
  --input "$CAPSACIN_INPUT" --work-dir /path/to/new-run \
  slice -- --symmetry 5 --weight 0.7 --auto
```

For unfamiliar inputs or ROI setup, use the inexpensive residue inspection to
check MODELs, chains, and actual residue ranges. It does not certify atom ordering
or assembly completeness. List candidate axes when the user asks to compare or
choose them, or when ranking needs inspection. Do not require an extra axis run
for every slice, or stop after inspection when slicing was requested.

The launcher copies the input into a new run directory, calls the original CLI,
streams and saves its output, and returns its exit status. It does not modify
the input or compute geometry itself. Supply a fresh `--work-dir` for each
axis/ROI/fold/weight run; without it, a unique `capsacin-*` directory is created
under the current directory and retained. Existing run directories are not
reused because native filenames omit ROI, axis rank, and manual reference values.

Read [references/cli.md](references/cli.md#selection-modes) for axis listing,
ROI/manual examples, parameter meanings, or optional restraint generation.

## Check and report the result

- Wait for process completion and read the exit status and `run.log`. On failure,
  report the actual cause and log path. Do not switch selection mode, change ROI,
  tune weight, or enable the legacy plane heuristic just to make the run pass.
- A successful slice normally writes
  `output/<input-stem>-sliced-sym<fold>-w<float-weight>.pdb` inside the run
  directory. An axis-only run writes a log and no sliced PDB. Inspect the actual
  output; the native pipeline can successfully write an empty structure.
- Report the input, mode, fold, weight, ROI/frame or reference indices as relevant,
  selected axis information available in the log, output path, and retained
  atom/residue/chain counts. The reported retained percentage uses the aligned
  pre-slice atom count as its denominator; it is not a mass percentage.
- Output coordinates are aligned/rebased and chains are renamed. Weight controls
  a geometric z cutoff; it is not a retained mass percentage. A sliced structure
  does not establish pore function or constitute a prepared MD system.
