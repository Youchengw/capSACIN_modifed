# capSACIN CLI reference

## Environment

Use the source repository at
<https://github.com/Youchengw/capSACIN_modifed>. The skill folder can be copied
elsewhere, but the engine and its data stay in a separate checkout.

For a new installation, the repository's setup is:

```bash
git clone https://github.com/Youchengw/capSACIN_modifed.git capSACIN
cd capSACIN/env
source create-env.sh
conda activate capSACIN
```

Run this in a shell with Conda initialized. `create-env.sh` creates the named
Python 3.9 environment, installs MDAnalysis from conda-forge, and installs the
runtime and development requirements. Reuse an existing environment when
available. Installing `requirements.txt` alone does not install MDAnalysis.
For a CLI-only environment, the development requirements are unnecessary.

Find the active interpreter with `python -c 'import sys; print(sys.executable)'`.
Use that path for `CAPSACIN_PYTHON` or `--python`; no desktop build is required.
For import or binary compatibility failures, fix the environment and report
that execution failed. Do not replace MDAnalysis calculations with mocked data.

## Selection modes

The examples below use the variables defined in `SKILL.md`. Values are examples,
not recommended scientific defaults. Each `slice` command creates a separate
retained run directory under the current directory. To choose the destination,
put `--work-dir /path/to/new-run` before `slice`.

### List axes only

```bash
"$CAPSACIN_PYTHON" "$CAPSACIN_RUNNER" --repo "$CAPSACIN_REPO" \
  --input "$CAPSACIN_INPUT" slice -- --symmetry 5 --list-axes
```

Use `--list-axes` without `--auto` to print candidates and exit. Combining them
runs slicing. Without an ROI, `--axis-index` refers to the global ranking.

### ROI-aware automatic slicing

Inspect the input's actual chain and residue ranges first. For the bundled
`1k3v.pdb`, `protein and chainid A and resid 250:290` is an example selection:

```bash
"$CAPSACIN_PYTHON" "$CAPSACIN_RUNNER" --repo "$CAPSACIN_REPO" \
  --input "$CAPSACIN_REPO/systemSetup/input/1k3v.pdb" inspect

"$CAPSACIN_PYTHON" "$CAPSACIN_RUNNER" --repo "$CAPSACIN_REPO" \
  --input "$CAPSACIN_REPO/systemSetup/input/1k3v.pdb" slice -- \
  --symmetry 5 --weight 0.7 --auto --axis-index 0 \
  --roi-selection "protein and chainid A and resid 250:290" --roi-frame 0
```

With an ROI, axis index 0 is the top **ROI-ranked** candidate. For inspection,
replace `--auto` with `--list-axes` while keeping the same selection and frame.
The ROI is evaluated in the specified MODEL, not across the whole assembly.
If no atoms match, correct the selection using actual residue IDs; do not
silently remove the ROI. In manual mode, `--roi-selection` is ignored.

### Manual reference atoms

The bundled `1k3v` fivefold manual example uses atom index 959:

```bash
"$CAPSACIN_PYTHON" "$CAPSACIN_RUNNER" --repo "$CAPSACIN_REPO" \
  --input "$CAPSACIN_REPO/systemSetup/input/1k3v.pdb" slice -- \
  --symmetry 5 --weight 0.7 --refindex 959
```

Use reference indices appropriate to the actual input. Add
`--legacy-plane-heuristic` only when the user requests that existing method.
Do not use it as a fallback after failure.

### Parameter and output details

| Native option | Meaning |
| --- | --- |
| `--symmetry`, `--sym` | 2, 3, or 5 |
| `--weight` | In [0, 1]; retains aligned atoms at `z >= max(z) * weight`, then applies existing cleanup |
| `--refindex` | One or more zero-based reference atom indices for manual mode |
| `--auto` | Automatic axis and reference selection |
| `--axis-index` | Zero-based rank in the candidate list for the current ROI/global mode |
| `--roi-selection` | Quoted MDAnalysis atom selection; `resid` is distinct from atom `index` |
| `--roi-frame` | Zero-based MODEL position |
| `--plot` | Matplotlib diagnostic plots; launcher defaults to headless `Agg` |

The launcher owns `--pdb` and derives its prefix from `--input`; do not pass a
second input prefix. It copies bytes, without converting formats or units.
It forwards other native options, including `--help`, unchanged. Wrapper options
must precede the command. To use an interactive plotting backend, explicitly set
`MPLBACKEND` to one supported by the environment.

The native CLI uses the invocation directory's `input/` and `output/`. If using
it directly, work from `systemSetup/` and pass a prefix such as `--pdb 1k3v`, not
an absolute path or a filename with `.pdb`. The launcher handles external input
paths by staging each run separately.

Logs and PDBs remain in the printed run directory, including after failure.
Weight in filenames is Python's float representation, e.g. `0.7` or `1.0`.
Different ROI/axis/reference runs can otherwise share a native output filename.
Read the reported paths instead of guessing a global `output/` location.

## Optional restraints

Only run `systemSetup/genRestraints.py` when restraints are requested. It is a
separate legacy CLI, not an MD simulation or topology builder. From a successful
slice's run directory, the command shape is:

```bash
cd /path/to/completed-run
MPLBACKEND=Agg "$CAPSACIN_PYTHON" "$CAPSACIN_REPO/systemSetup/genRestraints.py" \
  --pdb 1k3v --symmetry 5 --weight 0.7
```

Here `--weight` both chooses the sliced input filename and scales force constants
as `1000 * weight * s`. It is not independently adjustable from the input suffix
in the current CLI. The code hardcodes a full-weight region below z=15 Å and a
sigmoid region at 15–30 Å. State these semantics before treating the output as
the requested restraints. Different strength and slice-weight values need a
separate implementation change, not a silently renamed input.

Output is `output/<prefix>-w<weight>-posre/*chain_<id>.itp`; this directory omits
symmetry, so use the matching run directory. Indices restart at 1 for each chain
in the sliced PDB. Check them against the actual target topology before MD use.
