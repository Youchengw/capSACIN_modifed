# Auto ROI Candidate Discovery Plan

## Motivation

The current ROI-aware workflow improves on the legacy `--refindex` interface,
but it still requires the user to manually choose a residue range:

```bash
--roi-selection "protein and resid 100:140"
```

This is still tricky because:

- residue numbering differs for every capsid;
- a contiguous residue range is not always spatially compact;
- users may not know which residues lie near a 5-fold, 3-fold, or 2-fold axis;
- the original `refindex` workflow implicitly required expert knowledge of a
  useful atom in `MODEL 1`;
- `--roi-selection` makes the ROI wider and more explicit, but it does not
  automatically discover biologically meaningful capsid surface patches.

The next step should be to invert the current workflow:

```text
current: user ROI selection -> nearest symmetry axis
next:    symmetry axis -> nearby residue/atom ROI candidates
```

This lets the tool list candidate ROIs for every 5/3/2-fold axis, so the user
does not need to guess residue IDs.

## Current State

Implemented pieces:

- `findSymmetryAxes.build_icosahedral_axes(...)`
  - finds all candidate 5-fold, 3-fold, and 2-fold axes from MODEL COMs.
- `findSymmetryAxes.rank_axes_for_roi(...)`
  - ranks candidate axes by proximity to a user-supplied ROI center.
- `sliceCapsid.py --roi-selection ... --list-axes`
  - prints axis candidates ranked by a manually chosen ROI.
- `inspectResidueRanges.py`
  - reports valid residue ID ranges per PDB/MODEL/chain.

Main limitation:

- The user still needs to decide the ROI residue range before seeing useful
  axis-local residue suggestions.

## Proposed User Experience

Add a mode that lists axis-centered ROI candidates:

```bash
python sliceCapsid.py \
  --pdb 9jjh \
  --symmetry 5 \
  --list-axis-rois
```

Example output:

```text
[axis-roi] fold rank axis-index frame chain center-resid nearby-resids min-dist mean-dist axis
[axis-roi] 5    0    2          17    A     128          112:146      3.2      8.7       [...]
[axis-roi] 5    1    4          03    A     184          170:203      3.8      9.1       [...]
...
```

Then slicing can use either:

```bash
python sliceCapsid.py --pdb 9jjh --symmetry 5 --auto --axis-index 0
```

or a future explicit selector:

```bash
python sliceCapsid.py --pdb 9jjh --symmetry 5 --auto --axis-id 5fold:0
```

Optional convenience: print the generated ROI selection string:

```text
suggested-selection: "protein and name CA and resid 112:146"
suggested-roi-frame: 17
```

This bridges the new discovery mode with the already implemented
`--roi-selection` path.

## Core Algorithm

For each candidate symmetry axis:

1. Orient axis consistently.
   - For each physical axis, either sign is equivalent.
   - Choose the sign that points toward the nearest monomer/patch.

2. Compute residue-level representatives for all frames.
   - For each MODEL/frame and chain:
     - group protein atoms by `(frame, chainID, resid, resname)`;
     - compute residue center, preferably CA position when available;
     - compute radial distance from capsid center;
     - compute distance to the candidate axis line.

3. Rank residue representatives by distance to the axis line.
   - Line distance:

   ```text
   d = || cross(residue_center - capsid_center, axis) ||
   ```

   - Also track signed projection:

   ```text
   projection = dot(residue_center - capsid_center, axis)
   ```

   - Prefer positive projection after axis orientation so suggested ROI lies
     on the side that will become high-z after alignment.

4. Identify a compact surface patch.
   - Start from the nearest residue or nearest few residues.
   - Expand to residues on the same frame/chain that are:
     - spatially near the seed residue, and/or
     - contiguous in residue ID.
   - Suggested default:

   ```text
   take residues on the selected frame/chain with:
     distance_to_axis <= min_distance + 8 Å
     and radial_distance in the top 50% for that frame/chain
   ```

5. Convert selected residues into human-readable suggestions.
   - If mostly contiguous:

   ```text
   resid 112:146
   ```

   - If discontinuous:

   ```text
   resid 112 113 118 119 131 132
   ```

   - Prefer outputting both:
     - a compact range for convenience;
     - exact residue IDs for accuracy.

## Important Design Choices

### Axis-to-ROI is better than ROI-to-axis for discovery

The capsid has a known finite set of axes:

```text
5-fold: 6 axes
3-fold: 10 axes
2-fold: 15 axes
```

These are natural objects for automatic discovery. Once all axes are known,
the tool can ask: "what residue patch sits nearest each axis?"

### Residue IDs are only labels

Residue IDs should not drive the algorithm by themselves. They are useful for
display and for constructing MDAnalysis selection strings, but spatial
distance to an axis should drive discovery.

### Surface filtering matters

Without a surface/radial filter, the nearest residue to an axis may be buried
inside the protein. For CapSACIN surface slicing, candidates should bias toward
outer/surface residues.

Simple first implementation:

```text
surface_score = radial_distance / max_radial_distance_in_frame
keep residues with surface_score >= 0.5 or 0.6
```

This is crude but likely good enough for candidate listing.

Better future implementation:

- use solvent accessible surface area;
- use alpha shape / convex hull;
- use residue contact degree as a burial proxy.

### Frame choice should be explicit

Current `--roi-frame` defaults to `0`, but auto discovery should report the
frame it found:

```text
frame 17 means use --roi-frame 17 if routing through --roi-selection.
```

For direct `--axis-index` slicing, frame can be used internally when picking
reference atoms, but users should still see it in diagnostics.

## Proposed Functions

Add to `capsacin/findSymmetryAxes.py`:

```python
def compute_residue_representatives(universe, capsid_center):
    """
    Return a list of residue-level records:
      frame, chain_id, resid, resname, center, ca_center, radial_distance
    """

def rank_residues_for_axis(residue_records, axis, capsid_center,
                           prefer_positive_axis=True,
                           surface_quantile=0.5):
    """
    Return residue records ranked by distance to one axis.
    """

def summarize_axis_roi(axis, residue_records, capsid_center,
                       window_radius=8.0, max_residues=40):
    """
    Produce one suggested ROI patch for one axis.
    """

def list_axis_roi_candidates(universe, symmetry_type,
                             surface_quantile=0.5,
                             window_radius=8.0,
                             n_grid=2000):
    """
    Return one or more ROI suggestions for every candidate axis.
    """
```

Add to `sliceCapsid.py`:

```bash
--list-axis-rois
--surface-quantile 0.5
--roi-window-radius 8.0
```

Optional later:

```bash
--axis-id 5fold:0
--auto-roi
```

## Output Schema

Each axis ROI candidate should include:

```python
{
    "fold": 5,
    "axis_rank": 0,
    "axis_index": 2,
    "axis": np.ndarray(3),
    "symmetry_score": float,
    "frame": 17,
    "chain_id": "A",
    "center_resid": 128,
    "center_resname": "TYR",
    "min_axis_distance": float,
    "mean_axis_distance": float,
    "radial_distance": float,
    "residue_ids": [112, 113, ...],
    "residue_ranges": ["112:146"],
    "suggested_selection": "protein and name CA and resid 112:146",
    "suggested_roi_frame": 17,
}
```

## Validation Plan

### Smoke tests

Run:

```bash
python sliceCapsid.py --pdb 1k3v --symmetry 5 --list-axis-rois
python sliceCapsid.py --pdb 1k3v --symmetry 3 --list-axis-rois
python sliceCapsid.py --pdb 1k3v --symmetry 2 --list-axis-rois
python sliceCapsid.py --pdb 9jjh --symmetry 5 --list-axis-rois
```

Expected:

- 5-fold reports 6 candidates;
- 3-fold reports 10 candidates;
- 2-fold reports 15 candidates;
- each candidate reports a valid frame, chain, residue IDs, and selection.

### Compare against historical refindex examples

Known examples:

```bash
1k3v 5-fold refindex 959
1k3v 3-fold refindex 3131
1k3v 2-fold refindex 4073
```

Check:

1. map these atom indices to residues;
2. see whether those residues appear near one of the suggested axis ROI
   candidates;
3. confirm slicing with suggested ROI produces a surface similar to manual mode.

### VMD visual inspection

For each suggested ROI:

1. generate a full-capsid visualization PDB if needed;
2. color the suggested residues;
3. check whether the suggested patch lies near the intended symmetry feature.

## Risks and Caveats

- Nearest-to-axis residues may be buried, hence surface filtering is important.
- Residue ID ranges may not be spatially contiguous.
- Multi-chain MODELs, such as `4oq8` and `8des`, need chain-aware grouping.
- For 2-fold axes, the nearest patch may sit on an interface where two
  monomers are similarly close; output may need multiple patches per axis.
- If the PDB has missing loops, the nearest resolved residues may not represent
  the true biological ROI.

## Recommended Implementation Order

1. Add residue representative extraction.
2. Add axis-to-residue distance ranking.
3. Add `list_axis_roi_candidates()` returning structured dictionaries.
4. Add `sliceCapsid.py --list-axis-rois` printing a readable table.
5. Test on `1k3v`, `9jjh`, and one multi-chain MODEL input (`4oq8` or `8des`).
6. Add documentation to `roi_auto_axis_testing.md`.
7. Optional: wire suggested ROI directly into `--auto-roi`.

## Bottom Line

The current `--roi-selection` is useful but still manual.  The next useful
automation is not a smarter residue selection parser; it is an axis-centered ROI
discovery tool:

```text
Find all symmetry axes -> for each axis, suggest nearby surface residues.
```

That would let users test unfamiliar capsids without first guessing residue
numbers by hand.
