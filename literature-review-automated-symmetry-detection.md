# Automated Icosahedral Symmetry Axis Detection for CapSACIN

## A. Literature Review

### Paper 1: GRADE — "GRADE: A Code to Determine Clathrate Hydrate Structures"

**Citation:** Mahmoudinobar, F., Sarupria, S. et al. (2019). *Computer Physics Communications*, 244, 385–391.
DOI: `10.1016/j.cpc.2019.06.004` — PII: `S0010465519301833`

#### What It Does

GRADE identifies **clathrate hydrate cage structures** from atomic coordinates of water molecules (oxygen positions). It classifies water cages into types: 5¹² (pentagonal dodecahedron, 20 water molecules), 6²5¹² (tetrakaidecahedron, 24 waters), 6⁴5¹² (hexakaidecahedron, 28 waters), and 6⁸ (icosahedron-like, 36 waters).

The notation XⁿYᵐ means: n faces with X sides, m faces with Y sides. For example, 5¹² = 12 pentagonal faces, forming a dodecahedron.

#### Algorithm Detail

**Step 1 — Oxygen-oxygen connectivity graph:**

- Build an undirected graph where nodes are water oxygen atoms.
- Edge exists if O–O distance < 3.5 Å (first coordination shell cutoff).
- Each oxygen should have 4 neighbors (tetrahedral coordination in water).

**Step 2 — Ring detection (the core algorithm):**

- For each node, find all **pentagonal (5-membered) rings** and **hexagonal (6-membered) rings** it participates in.
- Ring-finding uses a **breadth-first search with path-memory**: starting from a node, traverse the connectivity graph tracking visited nodes; a ring is found when you return to a previously visited node at the correct step count (5 or 6).
- They apply a **filter** to only count *irreducible* rings (rings that are not composed of smaller rings).

**Step 3 — Cage classification:**

- For each oxygen, count how many 5-rings and 6-rings it belongs to.
- Match the ring count pattern to known clathrate cage signatures:
  - 5¹² cage: each oxygen belongs to exactly 3 five-rings, 0 six-rings
  - 6²5¹² cage: varied pattern (some O in 3 five-rings, some in 2 five-rings + 2 six-rings)
  - 6⁴5¹² cage: varied pattern
- Additional validation using the **F4 order parameter** (torsional angle distribution of H₂O pairs).

---

### Paper 2: CHILL+ — "Identification of Clathrate Hydrates, Hexagonal Ice, Cubic Ice, and Liquid Water in Simulations"

**Citation:** Nguyen, A. H. & Molinero, V. (2015). *J. Phys. Chem. B*, 119(29), 9369–9376.
DOI: `10.1021/jp510289t`

#### What It Does

CHILL+ classifies **individual water molecules** in MD simulations into one of: **cubic ice (Ic), hexagonal ice (Ih), clathrate hydrate, interfacial ice/clathrate, or liquid water**. It is an extension of the original CHILL algorithm (Moore et al., 2010) that adds clathrate detection.

#### Algorithm Detail

**Step 1 — Compute local bond order parameter q₃ for each water molecule:**

- For molecule *i*, identify its 4 nearest-neighbor O atoms (within 3.5 Å).
- For each pair of neighbors *(j, k)*, compute the coherence:
  ```
  c(i,j,k) = Σ_{m=-3}^{3} q̂_{3m}(i,j) · q̂_{3m}*(i,k)
  ```

  where q̂_{3m} is the normalized complex bond order parameter (Steinhardt order parameter with l=3), which captures the local orientational order of bonds around molecule *i*.
- c(i,j,k) ranges from −1 (perfectly eclipsed bond pair) to +1 (perfectly staggered bond pair).

**Step 2 — Classify each O–O bond as staggered or eclipsed:**

- **Staggered bond:** c(i,j,k)  =< −0.8 (bond pair is in staggered conformation, characteristic of cubic ice)
- **Eclipsed bond:** 0.25 >= c(i,j,k) >= 0.25 (bond pair is in eclipsed conformation, characteristic of clathrate — broader range than the original CHILL to capture the clathrate signal)
- The threshold values were empirically tuned to maximize detection rates near melting temperatures.

**Step 3 — Classify each water molecule by its bond-type profile:**


| # Staggered bonds | # Eclipsed bonds | Classification                                      |
| ----------------- | ---------------- | --------------------------------------------------- |
| 4                 | 0                | **Cubic ice (Ic)**                                  |
| 3                 | 1                | **Hexagonal ice (Ih)** — one eclipsed along c-axis |
| 2                 | 0                | **Interfacial ice**                                 |
| 0                 | 4                | **Clathrate hydrate**                               |
| any               | 3                | **Interfacial clathrate**                           |
| anything else     | anything else    | **Liquid water**                                    |

Additional criteria:

- Interfacial classifications require the molecule to have at least one neighbor classified as bulk ice or bulk clathrate (to suppress false positives in liquid).
- Liquid water is the default when no crystal signature matches.


---


## B. Current CapSACIN Workflow: How Reference Atoms Are Used

### The manual reference atom problem

The current workflow in `sliceCapsid.py` requires the user to provide `--refindex`, a VMD-format (0-based) atom index. This atom must be **in the Region of Interest (ROI)** near the desired symmetry axis. The algorithm then:

1. **Finds the reference atom (`pointA`):**

   ```python
   refPos = u.select_atoms(f"protein and index {indicesVMD[0]}")
   pointA = refPos.positions[0]
   ```
2. **Finds all symmetry-related copies** of that atom across the icosahedral capsid:

   ```python
   candidatePoints = chains[(chains.resids == refPos.resids[0]) & 
                            (chains.resname == refPos.resnames[0]) & 
                            (chains.name == refPos.names[0])]
   ```

   Since an icosahedral capsid has 60 asymmetric units (for T=1), there are ~60 copies of each atom. These copies encode the icosahedral symmetry.
3. **Selects `pointB` and `pointC` by nearest-neighbor distance** from `pointA`, with the neighbor rank depending on the desired fold type:

   ```
   2-fold: k1=0, k2=1  →  pointB = self (z+0.1), pointC = nearest copy
   3-fold: k1=1, k2=2  →  pointB = nearest copy, pointC = 2nd-nearest copy
   5-fold: k1=1, k2=3  →  pointB = nearest copy, pointC = 3rd-nearest copy
   ```

   A tiny z-perturbation (±0.1 Å) is added to break any degeneracy.
4. **Defines the symmetry axis normal** via cross product (in `definePlane.py`):

   ```python
   # u = vector from COM to pointB, v = vector from COM to pointC
   u_cross_v = np.cross(u, v)  # This is the symmetry axis direction
   ```
5. **Aligns the axis to the z-direction** using Rodrigues' rotation formula:

   ```python
   A = rotationMatrix(normalVector, [0, 0, 1])
   newCoords = coords @ A.T
   ```

### Why the distance-based heuristic works

For an icosahedral capsid, the ~60 symmetry-related copies of a given atom are arranged on the sphere. The nearest-neighbor distance rank correlates with the icosahedral symmetry group structure:

- For a **5-fold axis**: copies at indices k=1 and k=3 in the sorted distance list tend to be in adjacent asymmetric units related by 5-fold rotation.
- For a **3-fold axis**: copies at k=1 and k=2 are 3-fold-related.
- For a **2-fold axis**: copies at k=0 and k=1 are 2-fold-related.

This works **most of the time** but is fragile — it depends on the specific atom chosen, the PDB structure quality, and assumes the distance ordering reliably separates symmetry types.

### The real bottleneck

The manual step is not just entering an atom index — it requires the user to:

1. Open the capsid structure in VMD
2. Visually identify the region corresponding to a 5-fold, 3-fold, or 2-fold interface
3. Find a suitable atom in that region
4. Note its VMD index
5. Repeat for each symmetry type and each capsid

This is time-consuming and error-prone, especially for new capsid structures. The goal is to automate this entirely.

---

## C. Proposed Method: Automated Icosahedral Symmetry Axis Detection

### Core Idea

An icosahedral virus capsid (T=1) is composed of 60 identical protein chains (asymmetric units) arranged with perfect icosahedral (Iₕ) point group symmetry. The 60 chains are generated by applying the 60 rotation matrices of the icosahedral group to a single reference chain.

If we compute the **center of mass (COM) of each chain** and project them onto a unit sphere centered at the capsid COM, the resulting 60 points on the sphere form a pattern that **completely encodes** the icosahedral symmetry. The problem reduces to: **given 60 points on a sphere approximately obeying icosahedral symmetry, find the symmetry axes.**

The icosahedral group has:

- **6 five-fold axes** (through 12 vertices of an icosahedron): rotation by 72° leaves the structure invariant
- **10 three-fold axes** (through 20 face centers): rotation by 120°
- **15 two-fold axes** (through 30 edge midpoints): rotation by 180°

All three types of axes can be derived geometrically from the chain COM distribution.

---

### Method 1 (Primary — Inertia Tensor + Golden Ratio): Eigen-decomposition of the Chain COM Distribution

This is the recommended approach. It is simple, robust, and requires only elementary linear algebra.

#### Step 1.1: Compute the capsid center and chain COMs

```
Input: MDAnalysis Universe of the capsid
Output: capsid_center (3-vector), chain_coms (N×3 array)

1. capsid_center = mean position of all protein atoms
2. For each unique chain ID:
     chain_coms[i] = mean position of all atoms in that chain
   → Output: N chain COMs (N = 60 for T=1 capsid)
```

#### Step 1.2: Center the COMs and compute the 3×3 covariance matrix

```
1. centered_coms = chain_coms - capsid_center   # shape (N, 3)
2. Σ = (centered_coms.T @ centered_coms) / N    # 3×3 covariance matrix
```

#### Step 1.3: Eigendecomposition → three orthogonal 2-fold axes

```
eigenvalues, eigenvectors = np.linalg.eigh(Σ)
# eigenvectors are the principal axes, sorted by eigenvalue
# For a perfectly icosahedral distribution, all 3 eigenvalues are equal
# (the distribution is spherically symmetric), but the eigenvectors
# still provide a fixed orthogonal reference frame

# The 3 eigenvectors correspond to 3 mutually orthogonal 2-fold axes
# of the icosahedron. An icosahedron has 15 two-fold axes total;
# any 3 mutually orthogonal ones form a valid reference frame.
two_fold_axes = eigenvectors  # shape (3, 3), columns are axis directions
```

Why this works: The 60 chains of a T=1 capsid are arranged with icosahedral symmetry. The moment of inertia tensor (equivalently, the covariance matrix of chain COM positions) is isotropic for a perfect icosahedron (all eigenvalues equal), meaning the eigenvectors can be chosen as any orthonormal basis. However, for actual PDB structures with slight deviations from perfect symmetry, the eigenvectors will align with the directions of maximum variance — which correspond to 2-fold axes of the icosahedron, since these are the directions where chains appear in pairs.

#### Step 1.4: Build the icosahedral axis reference frame from 2-fold axes

Given 3 mutually orthogonal 2-fold axes {ê₁, ê₂, ê₃}, all symmetry axes of the icosahedron can be expressed as linear combinations. The standard construction uses the golden ratio φ = (1 + √5) / 2 ≈ 1.618.

The 6 five-fold axis directions (± each) are:

```
(±1, ±φ, 0)    in the {ê₁, ê₂, ê₃} basis (cyclic permutations)
(0, ±1, ±φ)
(±φ, 0, ±1)

This gives 12 vectors = 6 axes (each axis has two opposite directions).
```

The 10 three-fold axis directions (± each):

```
(±1, ±1, ±1)   in the {ê₁, ê₂, ê₃} basis
(±φ, ±1/φ, 0)  and cyclic permutations

This gives 20 vectors = 10 axes.
```

The 15 two-fold axis directions (± each):

```
ê₁, ê₂, ê₃ (the 3 principal axes)
(±φ, ±1, ±(φ-1)) / 2  and cyclic permutations

This gives 30 vectors = 15 axes.
```

All vectors are normalized to unit length.

#### Step 1.5: Select the desired symmetry axis

```
Given user-specified symmetry_type ∈ {2, 3, 5}:
  axes = compute_all_axes_of_type(symmetry_type, two_fold_axes)
  # For 5-fold: 6 candidate axes
  # For 3-fold: 10 candidate axes
  # For 2-fold: 15 candidate axes

Select one axis. Multiple strategies:
  a) Pick the axis with the most chains near it (densest cluster of chain COMs 
     near the axis direction on the sphere)
  b) Pick the axis whose positive hemisphere direction is closest to the 
     +z direction (for consistency)
  c) Return all axes and let the user choose by index (e.g., --axis-index 0)

Strategy (a) is recommended as default — it picks the most "canonical" 
axis of that type.
```

#### Step 1.6: Find a reference atom near the selected axis

```
1. For each chain, compute the angle between its COM direction and the axis:
     angle[i] = arccos(|chain_com_dir[i] · axis_dir|)
   
2. Find the chain with the smallest angle → this is the chain closest to the axis

3. Select an atom from this chain as the reference:
   - Use the Cα atom of a central residue (e.g., median residue by sequence)
   - Or use the chain's COM-proximal atom
   
4. Return this atom's index as the equivalent of --refindex
```

#### Step 1.7: Validate the output

```
1. Run the existing alignment algorithm with the auto-detected reference atom
2. Verify that the plane normal computed from the 3 reference points aligns
   with the symmetry axis to within a tolerance (~5°)
3. Optionally: compare alignment results with manually-curated reference indices
   for known capsids (1k3v)
```

---

### Method 2 (Alternative — Spherical K-Means / Density Peaks): Direct Axis Detection from Chain COM Distribution

If Method 1 proves insufficiently robust (e.g., for capsids with significant structural perturbations), this method finds symmetry axes directly by analyzing the spherical distribution of chain COMs.

#### Step 2.1: Project chain COMs onto the unit sphere

```
chain_dirs[i] = (chain_coms[i] - capsid_center) / ||chain_coms[i] - capsid_center||
```

#### Step 2.2: Find 5-fold axes by spherical density minima

On the sphere, a 5-fold axis passes through a point that is **maximally distant** from the nearest chain COM directions. This is because the 5-fold axis is at the center of a pentagonal arrangement of chains — the axis itself does not pass through any chain, but through the hole in the middle of 5 chains.

```
1. Generate a fine spherical grid (~10,000 points using Fibonacci sphere)
2. For each grid point, compute its angular distance to the nearest chain COM
3. Local maxima of "distance to nearest chain" are potential symmetry axes
4. Cluster these maxima → the 6 strongest clusters are the 5-fold axes
5. Verify: each 5-fold axis should have ~5 chains at similar angular distance (~37.4° 
   from the axis for a perfect icosahedron)
```

#### Step 2.3: Find 3-fold axes

A 3-fold axis passes through the center of a triangular face of the icosahedron. It is surrounded by 3 chains at equal angular distances (~41.8°).

```
1. For each grid point, find the 3 nearest chain COMs
2. Check if they are at approximately equal angular distances (within tolerance)
3. Check if the 3 chains are approximately 120° apart in azimuth around the axis
4. Cluster qualifying grid points → 10 clusters = 3-fold axes
```

#### Step 2.4: Find 2-fold axes

A 2-fold axis passes through the midpoint of an edge of the icosahedron, between 2 chains.

```
1. For each pair of nearby chains, compute the midpoint direction
2. Check if this direction is approximately equidistant from both chains
3. Verify that a 180° rotation around this axis maps the two chains to each other
4. Cluster qualifying directions → 15 clusters = 2-fold axes
```

#### Comparison: Method 1 vs. Method 2


| Criterion                                  | Method 1 (Inertia Tensor)                                | Method 2 (Spherical Density)                 |
| ------------------------------------------ | -------------------------------------------------------- | -------------------------------------------- |
| Complexity                                 | Low — eigendecomposition + analytic formulas            | Medium — spherical grid search + clustering |
| Robustness to noise                        | Good — eigendecomposition is stable                     | Good — density peaks are local and tolerant |
| Dependency on perfect icosahedral symmetry | Derives axes analytically from the group structure       | Finds axes empirically from data             |
| Computational cost                         | O(N) — trivial                                          | O(N_grid × N_chains) — small               |
| Handles incomplete capsids                 | Yes, if chain COM distribution approximately icosahedral | Yes, local density peaks still detectable    |
| Handles T>1 capsids                        | Yes (more chains, same symmetry)                         | Yes                                          |

**Recommendation:** Implement Method 1 first. It is simpler, mathematically clean, and leverages the known icosahedral group structure. Fall back to Method 2 only if Method 1 fails on specific capsids.

---

### Legacy Compatibility

The proposed implementation preserves full backward compatibility:

```
# Old (manual) usage — still works:
python sliceCapsid.py --pdb 1k3v --refindex 959 --symmetry 5 --weight 0.7

# New (auto) usage:
python sliceCapsid.py --pdb 1k3v --symmetry 5 --weight 0.7 --auto

# Auto with specific axis selection:
python sliceCapsid.py --pdb 1k3v --symmetry 5 --weight 0.7 --auto --axis-index 0
```

When `--auto` is passed, `--refindex` is ignored (or made optional). The code internally:

1. Calls the symmetry axis detection module
2. Finds the reference atom
3. Proceeds with the existing alignment + slicing pipeline unchanged

---

## D. Implementation Plan

### New file: `capsacin/findSymmetryAxes.py`

```
Functions:
├── compute_capsid_center(universe) → np.ndarray (3,)
│     Compute COM of all protein atoms
│
├── compute_chain_coms(universe) → np.ndarray (N×3)
│     Compute COM for each unique chain
│
├── find_principal_axes(chain_coms, capsid_center) → np.ndarray (3×3)
│     Eigendecomposition of covariance matrix → 3 orthogonal 2-fold axes
│
├── build_icosahedral_axes(two_fold_axes) → dict
│     Returns {'5fold': (6×3), '3fold': (10×3), '2fold': (15×3)}
│     Uses golden ratio formulas to derive all axes
│
├── select_best_axis(axes, chain_coms, capsid_center) → np.ndarray (3,)
│     Given a set of candidate axes, pick the one with best chain support
│
├── find_reference_atom(universe, axis, capsid_center) → int
│     Find an atom near the axis direction to use as reference
│
└── auto_detect_reference(pdb_path, symmetry_type) → int
      Top-level function: given a PDB and symmetry type, 
      return a reference atom index
```

### Modified file: `sliceCapsid.py`

- Add `--auto` flag and `--axis-index` to argument parser
- In `main()`: if `--auto`, call `auto_detect_reference()` instead of using `--refindex`
- The rest of the pipeline (pointB/pointC selection, plane normal, rotation, slicing) remains **unchanged**

### Verification Protocol

1. **Regression test on PPV (1k3v):**

   - Run auto-detection for 5-fold, 3-fold, and 2-fold
   - Compare resulting reference atom to manually curated values (959, 3131, 4073)
   - Verify the aligned structures match
2. **Cross-capsid validation:**

   - Test on all 13 capsids in `input/`
   - Verify each produces a valid sliced surface model
   - Check that symmetry-related chains are correctly identified
3. **Geometric validation:**

   - For each auto-detected axis, verify that rotating the capsid by 360°/n around that axis approximately reproduces the original structure
   - Quantify the RMSD between symmetry-related chains

---
