# ROI-Aware Automatic Axis Testing Guide

This guide is for testing the ROI-aware automatic symmetry-axis workflow on any
capsid in `systemSetup/input/`.

## 1. Environment Check

Run from the repository root:

```bash
cd /Users/yccccw/RESEARCH/capSACIN
conda activate capSACIN
cd systemSetup

python -c "import MDAnalysis; from capsacin import findSymmetryAxes; print('import ok')"
```

If this import fails, fix the `capSACIN` conda environment before testing
slicing. The base Python environment may have incompatible `numpy`/`h5py`
binaries.

## 2. Choose Test Parameters

For each capsid, fill in this table first:

| Variable | Meaning | Example |
|---|---|---|
| `PDB_ID` | Filename prefix in `input/{PDB_ID}.pdb` | `1k3v` |
| `FOLD` | Symmetry fold to test: `5`, `3`, or `2` | `5` |
| `ROI_SELECTION` | MDAnalysis selection for the ROI on one MODEL/frame | `"protein and resid 250:290"` |
| `ROI_FRAME` | MODEL/frame containing that ROI copy | `0` |
| `WEIGHT` | CapSACIN slicing weight | `0.63` |

First inspect the available residue IDs:

```bash
python inspectResidueRanges.py PDB_ID
```

Example:

```bash
python inspectResidueRanges.py 9jjh
```

Example output:

```text
PDB: /Users/yccccw/RESEARCH/capSACIN/systemSetup/input/9jjh.pdb
Models: 60

Global residue ranges by chain:
chain  min_resid  max_resid  unique_resids  insertion_codes
    A         66        233            168  -

Per-MODEL residue ranges:
model  roi_frame  chain  min_resid  max_resid  unique_resids
    1          0      A         66        233            168
```

This means valid ROI selections for `9jjh` should use residue IDs between
`66` and `233`, for example:

```bash
--roi-selection "protein and resid 100:140" --roi-frame 0
```

If you choose a range outside the listed residue IDs, for example
`"protein and resid 900:1100"` for `9jjh`, the ROI will match no atoms.

Good first-pass ROI selections:

```bash
"protein and resid 250:290"
"protein and resid 300:340"
"protein and name CA and resid 250:290"
```

Prefer a compact ROI near the expected symmetry feature. If the ROI is far from
all candidate axes, the tool still chooses the closest geometric match, but the
result is less meaningful and should be inspected with `--list-axes`.

## 3. Confirm the ROI Selection Matches Atoms

Use `--list-axes` first. This tests PDB loading, ROI parsing, axis detection,
and ROI ranking without writing a sliced PDB.

```bash
python sliceCapsid.py \
  --pdb PDB_ID \
  --symmetry FOLD \
  --roi-selection "ROI_SELECTION" \
  --roi-frame ROI_FRAME \
  --list-axes
```

Example:

```bash
python sliceCapsid.py \
  --pdb 1k3v \
  --symmetry 5 \
  --roi-selection "protein and resid 250:290" \
  --roi-frame 0 \
  --list-axes
```

Expected output format:

```text
[axes] rank axis-index roi-score roi-angle-deg roi-line-dist symmetry-score axis
[axes]    0 ...
[axes]    1 ...
```

Interpretation:

| Column | Meaning |
|---|---|
| `rank` | ROI-ranked candidate order; `0` is the default auto choice |
| `axis-index` | Original symmetry-ranked axis index |
| `roi-score` | Combined ROI proximity score; lower is better |
| `roi-angle-deg` | Angle between ROI direction and candidate axis |
| `roi-line-dist` | Distance from ROI center to the candidate axis line |
| `symmetry-score` | Global rotational symmetry quality; lower is better |

If `roi-angle-deg` and `roi-line-dist` are both large for every axis, try a
more localized or more appropriate ROI.

## 4. Run ROI-Aware Auto Slicing

Use a nonstandard `WEIGHT` during testing so you do not overwrite existing
reference outputs.

```bash
python sliceCapsid.py \
  --pdb PDB_ID \
  --symmetry FOLD \
  --weight WEIGHT \
  --auto \
  --roi-selection "ROI_SELECTION" \
  --roi-frame ROI_FRAME
```

Example:

```bash
python sliceCapsid.py \
  --pdb 1k3v \
  --symmetry 5 \
  --weight 0.63 \
  --auto \
  --roi-selection "protein and resid 250:290" \
  --roi-frame 0
```

Expected diagnostic lines:

```text
[auto] ROI selection: ...
[auto] Using ROI-ranked 5-fold axis at rank 0
[auto] Detected 5-fold axis: [...]
[auto] Reference atom indices: [...]
[auto] Reference frames: [...]
[auto] Using detected axis as alignment normal (bypassing 3-point method)
```

Expected output file:

```bash
ls -lh output/PDB_ID-sliced-wWEIGHT.pdb
```

Example:

```bash
ls -lh output/1k3v-sliced-w0.63.pdb
```

## 5. Test Another Candidate Axis

After reviewing `--list-axes`, choose a different ROI rank with `--axis-index`.
In ROI-aware mode, `--axis-index` refers to the printed `rank`, not the original
`axis-index` column.

```bash
python sliceCapsid.py \
  --pdb PDB_ID \
  --symmetry FOLD \
  --weight 0.631 \
  --auto \
  --roi-selection "ROI_SELECTION" \
  --roi-frame ROI_FRAME \
  --axis-index 1
```

Compare the output with the rank-0 model to confirm that axis selection changes
the retained surface.

## 6. Test All Fold Types

For each capsid/ROI, run:

```bash
python sliceCapsid.py --pdb PDB_ID --symmetry 5 --roi-selection "ROI_SELECTION" --roi-frame ROI_FRAME --list-axes
python sliceCapsid.py --pdb PDB_ID --symmetry 3 --roi-selection "ROI_SELECTION" --roi-frame ROI_FRAME --list-axes
python sliceCapsid.py --pdb PDB_ID --symmetry 2 --roi-selection "ROI_SELECTION" --roi-frame ROI_FRAME --list-axes
```

Then slice the fold types you care about:

```bash
python sliceCapsid.py --pdb PDB_ID --symmetry 5 --weight 0.635 --auto --roi-selection "ROI_SELECTION" --roi-frame ROI_FRAME
python sliceCapsid.py --pdb PDB_ID --symmetry 3 --weight 0.636 --auto --roi-selection "ROI_SELECTION" --roi-frame ROI_FRAME
python sliceCapsid.py --pdb PDB_ID --symmetry 2 --weight 0.637 --auto --roi-selection "ROI_SELECTION" --roi-frame ROI_FRAME
```

Note: 2-fold auto mode preserves the legacy CapSACIN convention where `pointB`
comes from the same monomer as `pointA` and is separated by a small z
perturbation before alignment.

## 7. Regression Test Legacy Modes

The original manual reference-atom workflow is still supported. Test it with a
known example:

```bash
python sliceCapsid.py \
  --pdb 1k3v \
  --symmetry 5 \
  --refindex 959 \
  --weight 0.638
```

Also test global auto mode without ROI:

```bash
python sliceCapsid.py \
  --pdb 1k3v \
  --symmetry 5 \
  --weight 0.639 \
  --auto
```

These confirm that ROI-aware auto did not break the existing manual or global
auto paths.

## 8. Suggested Capsid Test Matrix

Start with small smoke tests on a few structures:

| PDB_ID | Fold | ROI selection | First command |
|---|---:|---|---|
| `1k3v` | 5 | `"protein and resid 250:290"` | `--list-axes` |
| `1k3v` | 3 | `"protein and resid 250:290"` | `--list-axes` |
| `1k3v` | 2 | `"protein and resid 250:290"` | `--list-axes` |
| `9jjh` | 5 | Choose a compact residue range from frame 0 | `--list-axes` |
| `3ra2` | 5 | Choose a compact residue range from frame 0 | `--list-axes` |

For new capsids, first inspect the residue numbering in `input/{PDB_ID}.pdb`,
then choose a compact residue window that exists in MODEL/frame 0.

## 9. Common Failures

| Symptom | Likely cause | Fix |
|---|---|---|
| `ROI selection matched no atoms` | Residue range or selection does not exist on `ROI_FRAME` | Adjust `ROI_SELECTION` or `ROI_FRAME` |
| `Found only ... monomers` | PDB is not a multi-MODEL capsid in the expected format | Check input formatting |
| No output PDB | Slicing/cleanup removed everything or save failed | Try lower `--weight`; inspect terminal messages |
| Axis ranking looks ambiguous | ROI is broad or far from all axes | Use a smaller ROI or inspect another fold |

## 10. Minimal Copy/Paste Template

```bash
cd /Users/yccccw/RESEARCH/capSACIN/systemSetup

PDB_ID=1k3v
FOLD=5
ROI='protein and resid 250:290'
ROI_FRAME=0
WEIGHT=0.63

python sliceCapsid.py --pdb "$PDB_ID" --symmetry "$FOLD" \
  --roi-selection "$ROI" --roi-frame "$ROI_FRAME" --list-axes

python sliceCapsid.py --pdb "$PDB_ID" --symmetry "$FOLD" \
  --weight "$WEIGHT" --auto \
  --roi-selection "$ROI" --roi-frame "$ROI_FRAME"

ls -lh "output/${PDB_ID}-sliced-w${WEIGHT}.pdb"
```
