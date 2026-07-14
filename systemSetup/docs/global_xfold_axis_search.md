# Global X-Fold Symmetry Axis Search

This document explains the global x-fold axis search implemented in
`systemSetup/capsacin/findSymmetryAxes.py`.

Here, "x-fold" means a rotational symmetry axis of order `x`, where the
current CapSACIN workflow supports:

- `5-fold`: 6 axes in an icosahedral capsid
- `3-fold`: 10 axes in an icosahedral capsid
- `2-fold`: 15 axes in an icosahedral capsid

The key idea is simple: instead of asking the user to manually pick a
reference atom, the code searches the whole sphere for directions that behave
like real rotational symmetry axes.

## Why This Code Exists

The original CapSACIN workflow needs a reference atom index, usually selected
manually in VMD. That manual choice encodes two things:

1. Which physical symmetry axis should be used.
2. Which side of the capsid should be kept after alignment and slicing.

The automatic workflow separates these concerns. The global axis search first
detects all candidate `5-fold`, `3-fold`, and `2-fold` axes from the whole
capsid geometry. Later code can either choose the best global axis or rank axes
by proximity to a user-defined ROI.

## Input Assumption

CapSACIN input PDB files are multi-MODEL capsid files. In MDAnalysis, those
`MODEL` records appear as trajectory frames:

```text
MODEL 1 -> universe.trajectory[0] -> one asymmetric unit copy
MODEL 2 -> universe.trajectory[1] -> another symmetry copy
...
MODEL 60 -> universe.trajectory[59]
```

The global axis search therefore treats each frame as one monomer or asymmetric
unit copy. It computes one center per frame, then uses those 60 centers as a
geometric proxy for the icosahedral capsid.

## High-Level Pipeline

The core pipeline is:

```text
MDAnalysis Universe
  -> compute capsid center
  -> compute one monomer center per MODEL/frame
  -> convert monomer centers into unit directions from the capsid center
  -> sample candidate axis directions on the unit sphere
  -> score each candidate by rotational self-consistency
  -> keep local minima in the score landscape
  -> remove duplicate +/- versions of the same physical axis
  -> return the best expected number of axes for the requested fold
```

The main functions are:

```python
compute_capsid_center(universe)
compute_monomer_coms(universe)
find_symmetry_axes(monomer_coms, capsid_center, fold, n_grid, n_expected)
build_icosahedral_axes(monomer_coms, capsid_center, n_grid)
```

## 1. Computing the Capsid Center

Key code:

```python
def compute_capsid_center(universe):
    all_positions = []
    for ts in universe.trajectory:
        protein = universe.select_atoms("protein")
        all_positions.append(protein.positions)
    stacked = np.vstack(all_positions)
    return stacked.mean(axis=0)
```

Purpose:

This computes one global geometric center for the entire capsid. Because each
trajectory frame is one capsid copy, the function iterates over all frames and
stacks all protein coordinates before averaging them.

Method:

- `universe.trajectory` loops over every MODEL/frame.
- `select_atoms("protein")` extracts protein atoms from the current copy.
- `np.vstack(all_positions)` creates one large coordinate array.
- `stacked.mean(axis=0)` returns the global center as an `(x, y, z)` vector.

Why it matters:

All later axis calculations need directions from the capsid center to each
monomer center. If the center is wrong, the symmetry score becomes noisy or
biased toward an off-center direction.

## 2. Computing Monomer Centers

Key code:

```python
def compute_monomer_coms(universe):
    coms = []
    for ts in universe.trajectory:
        protein = universe.select_atoms("protein")
        coms.append(protein.positions.mean(axis=0))
    return np.array(coms)
```

Purpose:

This returns one representative point for each capsid copy.

Method:

- Each trajectory frame is treated as one monomer/asymmetric unit.
- The function averages protein atom coordinates within that frame.
- The result is usually a `(60, 3)` array for a T=1 icosahedral capsid with 60
  symmetry copies.

Why it matters:

The global axis search does not use every atom directly. It first reduces the
capsid into 60 monomer center directions. That makes the symmetry search fast
and robust to local atomic detail.

## 3. Normalizing Monomer Directions

Key code from `find_symmetry_axes()`:

```python
monomer_dirs = monomer_coms - capsid_center
monomer_dirs = np.array([_normalize(d) for d in monomer_dirs])
```

Purpose:

This projects monomer centers onto a unit sphere centered at the capsid center.

Method:

- Subtracting `capsid_center` converts coordinates into radial vectors.
- `_normalize(d)` keeps only direction, not distance.

Why it matters:

Rotational symmetry is about angular relationships around the capsid. The exact
radial distance of a monomer center is less important than the direction in
which that monomer lies.

## 4. Sampling Candidate Axis Directions

Key code:

```python
def _fibonacci_sphere(n_samples):
    if n_samples < 2:
        raise ValueError(f"n_samples must be >= 2, got {n_samples}.")
    golden = (1.0 + np.sqrt(5.0)) / 2.0
    dirs = np.empty((n_samples, 3))
    for i in range(n_samples):
        y = 1.0 - (i / float(n_samples - 1)) * 2.0
        radius = np.sqrt(1.0 - y * y)
        theta = 2.0 * np.pi * i / golden
        dirs[i, 0] = np.cos(theta) * radius
        dirs[i, 1] = y
        dirs[i, 2] = np.sin(theta) * radius
    return dirs
```

Purpose:

This generates approximately uniform candidate directions on the unit sphere.

Method:

- `y` moves from the north pole to the south pole.
- `radius = sqrt(1 - y*y)` gives the circle radius at that height.
- `theta = 2*pi*i/golden` advances around the sphere using the golden ratio.
- Each candidate is stored as a unit vector.

Why Fibonacci sampling:

A latitude-longitude grid clusters points near the poles. Random sampling is
not reproducible. A Fibonacci sphere gives deterministic, reasonably even
coverage without extra dependencies.

Important parameter:

`n_grid` controls the number of candidate directions. Higher values improve
angular resolution but increase runtime linearly.

## 5. Scoring One Candidate Axis

The most important function is `_rotational_symmetry_score()`.

Key code:

```python
def _rotational_symmetry_score(axis, monomer_dirs, fold):
    axis = _normalize(axis)
    N = len(monomer_dirs)
    total_err = 0.0

    for k in range(1, fold):
        angle = k * 2.0 * np.pi / fold
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)

        dot_products = monomer_dirs @ axis
        cross_products = np.cross(axis, monomer_dirs)
        rotated = (
            cos_a * monomer_dirs
            + sin_a * cross_products
            + (1.0 - cos_a) * np.outer(dot_products, axis)
        )
```

Purpose:

This asks: if this candidate direction were a true `fold`-fold axis, would
rotating every monomer direction around it map the capsid back onto itself?

Method:

- For a `5-fold` axis, the code tests rotations by 72, 144, 216, and 288
  degrees.
- For a `3-fold` axis, it tests rotations by 120 and 240 degrees.
- For a `2-fold` axis, it tests one 180 degree rotation.
- Rodrigues' rotation formula is applied to all monomer direction vectors at
  once.

The vectorized Rodrigues formula is:

```text
rotated = cos(a) * v
        + sin(a) * (axis cross v)
        + (1 - cos(a)) * dot(v, axis) * axis
```

In the implementation:

- `dot_products = monomer_dirs @ axis` computes `dot(v, axis)` for every
  monomer.
- `np.cross(axis, monomer_dirs)` computes `axis cross v` for every monomer.
- `np.outer(dot_products, axis)` rebuilds the axis-parallel component.

Why it matters:

This directly tests rotational symmetry. No assumption is made about how the
capsid is oriented in the input PDB.

## 6. Matching Rotated Monomers Back to Real Monomers

Key code:

```python
norms = np.linalg.norm(rotated, axis=1, keepdims=True)
rotated = rotated / norms

cos_dists = np.abs(rotated @ monomer_dirs.T)
cos_dists = np.clip(cos_dists, -1.0, 1.0)
min_angles = np.min(np.arccos(cos_dists), axis=1)
total_err += np.sum(min_angles)
```

Purpose:

After rotating all monomer directions, the function finds the closest real
monomer direction for each rotated direction.

Method:

- `rotated @ monomer_dirs.T` creates an `(N, N)` cosine similarity matrix.
- Entry `(i, j)` measures how closely rotated monomer `i` matches original
  monomer `j`.
- `np.abs(...)` treats `axis` and `-axis` as the same physical line.
- `arccos` converts cosine similarity into angular error.
- For each rotated monomer, the smallest angular error is used.

Why lower score is better:

The returned value is the mean angular mismatch:

```python
return total_err / (N * (fold - 1))
```

A perfect symmetry axis would have score near `0.0`. Larger values mean the
candidate axis does not rotate the monomer distribution cleanly onto itself.

## 7. Finding Local Minima on the Sphere

Key code:

```python
def _find_local_minima(scores, directions, min_angle=np.deg2rad(5.0)):
    cos_thresh = np.cos(min_angle)
    minima = []
    for i in range(len(scores)):
        cos_dists = np.abs(directions @ directions[i])
        neighbours = np.where(cos_dists > cos_thresh)[0]
        if all(scores[i] <= scores[j] for j in neighbours):
            minima.append((scores[i], directions[i].copy()))
    minima.sort(key=lambda x: x[0])
    return minima
```

Purpose:

This turns a list of raw scores into candidate axis peaks. Since lower score is
better, a good axis is a local minimum in the spherical score landscape.

Method:

- For every sampled direction, find nearby sampled directions within
  `min_angle`.
- If the current point has a score no worse than all nearby points, keep it as
  a local minimum.
- Sort all minima by score.

Why `np.abs(directions @ directions[i])` is used:

A physical axis has no inherent sign. The directions `axis` and `-axis`
describe the same line, so the code compares directions using the absolute dot
product.

## 8. Deduplicating Physical Axes

Key code:

```python
def _deduplicate_minima(minima, min_separation=np.deg2rad(10.0)):
    cos_sep = np.cos(min_separation)
    unique = []
    for score, direction in minima:
        is_dup = False
        for _, existing in unique:
            if abs(np.dot(direction, existing)) > cos_sep:
                is_dup = True
                break
        if not is_dup:
            unique.append((score, direction))
    return unique
```

Purpose:

This removes repeated detections of the same physical axis.

Method:

- The minima are already sorted best-first.
- For each new minimum, compare it with axes already kept.
- If it is nearly collinear with an existing axis, skip it.
- Otherwise, keep it.

Why absolute dot product:

`axis` and `-axis` are the same physical axis. `abs(dot(...))` makes both signs
collapse into one candidate.

Why this is needed:

The spherical grid can place several sampled directions near the same true
axis. Local minima detection reduces this, but deduplication is the final guard
against returning the same axis multiple times.

## 9. Searching One Fold Type

Key code:

```python
def find_symmetry_axes(monomer_coms, capsid_center, fold,
                       n_grid=5000, n_expected=None):
    monomer_dirs = monomer_coms - capsid_center
    monomer_dirs = np.array([_normalize(d) for d in monomer_dirs])

    directions = _fibonacci_sphere(n_grid)

    scores = np.array([
        _rotational_symmetry_score(d, monomer_dirs, fold)
        for d in directions
    ])

    minima = _find_local_minima(scores, directions)
    unique = _deduplicate_minima(minima)

    if n_expected is not None:
        unique = unique[:n_expected]

    axes = np.array([d for _, d in unique])
    return axes
```

Purpose:

This is the main global x-fold search function. It finds axes for one requested
fold type.

Method:

1. Convert monomer centers to unit directions.
2. Generate candidate axis directions.
3. Score every candidate direction.
4. Keep local score minima.
5. Deduplicate equivalent physical axes.
6. Optionally keep only the expected number of axes.

Expected axis counts:

```text
5-fold -> 6 axes
3-fold -> 10 axes
2-fold -> 15 axes
```

Output:

The return value is an array of axis directions with shape `(K, 3)`. The axes
are sorted by increasing symmetry score, so index `0` is the best global match.

## 10. Building All Icosahedral Axes

Key code:

```python
def build_icosahedral_axes(monomer_coms, capsid_center, n_grid=2000):
    fivefold = find_symmetry_axes(
        monomer_coms, capsid_center, fold=5,
        n_grid=n_grid, n_expected=6,
    )
    threefold = find_symmetry_axes(
        monomer_coms, capsid_center, fold=3,
        n_grid=n_grid, n_expected=10,
    )
    twofold = find_symmetry_axes(
        monomer_coms, capsid_center, fold=2,
        n_grid=n_grid, n_expected=15,
    )
    return {"5fold": fivefold, "3fold": threefold, "2fold": twofold}
```

Purpose:

This convenience function runs the global search for all supported fold types.

Method:

- It calls `find_symmetry_axes()` three times.
- It uses the known icosahedral axis counts as `n_expected`.
- It returns a dictionary keyed by fold type.

Why `n_grid=2000` here:

The all-axis search is called inside the automatic pipeline, so the default is
set lower than the standalone `find_symmetry_axes()` default. This keeps the
workflow fast while still giving stable candidate axes for typical capsid
inputs.

## 11. Computing Search Context Once

Key code:

```python
def _compute_axis_context(universe, n_grid=2000):
    capsid_center = compute_capsid_center(universe)
    monomer_coms = compute_monomer_coms(universe)

    if len(monomer_coms) < 3:
        raise ValueError(...)

    monomer_dirs = monomer_coms - capsid_center
    monomer_dirs = np.array([_normalize(d) for d in monomer_dirs])
    all_axes = build_icosahedral_axes(monomer_coms, capsid_center, n_grid=n_grid)

    return capsid_center, monomer_coms, monomer_dirs, all_axes
```

Purpose:

This helper bundles the expensive shared setup used by automatic reference
detection, axis listing, and ROI-aware ranking.

Method:

- Compute the capsid center.
- Compute monomer centers.
- Validate that there are enough frames to infer symmetry.
- Compute normalized monomer directions.
- Build all 5-fold, 3-fold, and 2-fold axes.

Why it matters:

Without this helper, each caller would repeat the same setup and risk drifting
into slightly different behavior.

## 12. ROI-Aware Axis Ranking

The global search finds axes based only on symmetry quality. ROI ranking is a
separate layer on top of those global candidates.

Key code:

```python
roi_vector = roi_center - capsid_center
roi_dir = roi_vector / roi_distance

for axis_index, axis in enumerate(candidates):
    axis = _normalize(axis)
    if np.dot(axis, roi_dir) < 0:
        axis = -axis

    axis_dot = np.clip(np.dot(axis, roi_dir), -1.0, 1.0)
    roi_angle = np.arccos(axis_dot)
    roi_line_distance = np.linalg.norm(np.cross(roi_vector, axis))
    roi_distance_score = roi_line_distance / capsid_radius
    symmetry_score = _rotational_symmetry_score(axis, monomer_dirs, symmetry_type)

    roi_score = roi_angle + roi_distance_score + 0.25 * symmetry_score
```

Purpose:

This ranks already-detected global axes by how close they are to a user-defined
region of interest.

Method:

- Compute the direction from capsid center to ROI center.
- Flip each candidate axis so its positive direction points toward the ROI.
- Measure angular distance between the ROI direction and the axis.
- Measure line distance from the ROI center to the candidate axis.
- Add a small symmetry-score term so poor global axes are not preferred during
  ties.

Important behavior:

Without ROI selection, `axis_index=0` means the best global symmetry-ranked
axis. With ROI selection, `axis_index=0` means the closest ROI-ranked axis.

## 13. Automatic Reference Detection

The global x-fold search feeds into `auto_detect_reference()`.

Key code:

```python
capsid_center, monomer_coms, monomer_dirs, all_axes = _compute_axis_context(
    universe
)

key = _axis_key(symmetry_type)
candidates = all_axes[key]

if roi_selection:
    ranked_axes = _rank_axes_for_roi_from_context(...)
    best_axis = ranked_axes[axis_index]["axis"]
else:
    best_axis = candidates[axis_index]
```

Purpose:

This chooses the final axis used by `sliceCapsid.py --auto`.

Method:

- Compute all global axes.
- Select the requested fold type.
- If no ROI is provided, select by global symmetry rank.
- If ROI is provided, select by ROI rank.

After choosing `best_axis`, the function finds reference frames:

```python
ref_frames = _find_reference_frames(
    monomer_dirs,
    best_axis,
    symmetry_type,
    prefer_positive_axis=bool(roi_selection),
)
```

Then it chooses representative atoms from those frames:

```python
ref_indices = _pick_reference_atoms(
    universe, ref_frames, roi_selection=roi_selection
)
```

Output:

```python
return ref_indices, best_axis, ref_frames
```

These values tell `sliceCapsid.py` which atom coordinates to use for local
diagnostics and which global axis to use for alignment.

## 14. How `sliceCapsid.py` Uses the Global Axis

Key code from `sliceCapsid.py`:

```python
if auto_mode:
    auto_ref_indices, axis_dir, auto_frames = findSymmetryAxes.auto_detect_reference(
        u,
        symmetry,
        axis_index=axis_index,
        roi_selection=roi_selection,
        roi_frame=roi_frame,
    )
```

Purpose:

This connects the global axis search to the command-line `--auto` workflow.

Later, automatic mode uses the detected axis directly:

```python
if auto_mode and auto_ref_indices is not None:
    normalVector = axis_dir
    print("[auto] Using detected axis as alignment normal (bypassing 3-point method)")
```

Why this is important:

The old manual workflow estimates a plane normal from three selected points.
The automatic workflow uses the globally detected symmetry axis as the
alignment normal. That is usually more stable because it comes from all monomer
centers rather than a small local atom set.

## Practical Usage

List global candidate axes:

```bash
python systemSetup/sliceCapsid.py --pdb 1k3v --symmetry 5 --list-axes
```

Use the best global 5-fold axis:

```bash
python systemSetup/sliceCapsid.py --pdb 1k3v --symmetry 5 --auto
```

List ROI-ranked 5-fold axes:

```bash
python systemSetup/sliceCapsid.py \
  --pdb 1k3v \
  --symmetry 5 \
  --roi-selection "protein and resid 250:290" \
  --list-axes
```

Slice using the top ROI-ranked 5-fold axis:

```bash
python systemSetup/sliceCapsid.py \
  --pdb 1k3v \
  --symmetry 5 \
  --auto \
  --roi-selection "protein and resid 250:290"
```

## Debugging Notes

Useful outputs:

- `score`: global rotational symmetry score. Lower is better.
- `axis_index`: index in the global symmetry-ranked candidate list.
- `roi_rank`: index in the ROI-ranked list, only present when ROI selection is
  used.
- `roi_angle_deg`: angle between the ROI direction and the axis direction.
- `roi_line_distance`: shortest distance from ROI center to the axis line.

Common confusion:

- `resid` is a residue label used in MDAnalysis selections.
- `index` is an atom index.
- `MODEL` number is one-based in the PDB file, but `roi_frame` is zero-based in
  MDAnalysis. Therefore `MODEL 1` corresponds to `roi_frame=0`.

## Summary

The global x-fold search is a direct spherical search:

1. Reduce the multi-MODEL capsid to one center point per symmetry copy.
2. Project those centers onto a sphere.
3. Sample possible axis directions with a Fibonacci grid.
4. Rotate all monomer directions around each candidate axis.
5. Score how well the rotated directions match the original directions.
6. Keep local score minima and remove duplicate axes.
7. Return the expected number of 5-fold, 3-fold, or 2-fold axes.

This method is orientation-independent, works uniformly for 2-fold, 3-fold, and
5-fold symmetry, and gives a continuous quality score that can be used for
diagnostics and ranking.
