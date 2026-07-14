# CapSACIN — Capsid Surface Abstraction and Computationally-Induced Nanofragmentation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/)
[![Release v0.1.0](https://img.shields.io/badge/release-v0.1.0-2ea44f.svg)](https://github.com/Youchengw/capSACIN_modifed/releases/tag/v0.1.0)
[![macOS Apple Silicon](https://img.shields.io/badge/macOS-Apple%20Silicon-black.svg)](https://github.com/Youchengw/capSACIN_modifed/releases/download/v0.1.0/capSACIN-Studio-v0.1.0-macOS-arm64.dmg)

A computational framework for constructing atomistic surface models of icosahedral virus capsids, enabling high-throughput molecular dynamics simulations of virus–excipient interactions without the prohibitive cost of simulating fully assembled capsids.

---

## capSACIN Studio

**capSACIN Studio** is the graphical desktop interface for the CapSACIN workflow. It combines symmetry-axis selection and surface slicing controls with an integrated [Mol*](https://molstar.org/) 3D capsid viewer.

### Download

The current release is **v0.1.0** for Apple Silicon Macs running macOS 13 or later.

- [Download the macOS DMG](https://github.com/Youchengw/capSACIN_modifed/releases/download/v0.1.0/capSACIN-Studio-v0.1.0-macOS-arm64.dmg) — recommended installer.
- [Download the zipped application](https://github.com/Youchengw/capSACIN_modifed/releases/download/v0.1.0/capSACIN-Studio-v0.1.0-macOS-arm64.zip).
- [View the v0.1.0 release notes and checksums](https://github.com/Youchengw/capSACIN_modifed/releases/tag/v0.1.0).

The release bundles the Python sidecar and 13 example PDB structures. A separate Python or Conda environment is not required to run the packaged application.

### Install on macOS

1. Download and open the `.dmg` file.
2. Drag **capSACIN Studio** into `Applications`.
3. Launch the app and select a built-in capsid, or open a local PDB file.

The v0.1.0 build is ad-hoc signed and is not notarized with an Apple Developer ID. If macOS blocks the first launch, right-click the app and choose **Open**, or allow it from **System Settings → Privacy & Security**.

### Desktop workflow

1. Select a built-in structure or open a local PDB file.
2. Choose a 2-fold, 3-fold, or 5-fold symmetry axis.
3. Adjust the slicing weight `ω`; larger values remove more of the capsid.
4. Click **Prepare Preview** to detect and rank candidate symmetry axes.
5. Inspect the selected axis, slicing plane, and ROI in the 3D viewer.
6. Click **Run capSACIN** to compute the sliced structure.
7. Switch between **Original**, **Sliced**, and **Overlay** views, then save the resulting PDB.

Advanced settings expose the axis candidate rank, ROI chain and residue range, raw MDAnalysis selection overrides, and manual reference-index workflow. The legacy command-line path remains available and is documented below.

---

## Overview

Viral biologics — vaccines, gene therapy vectors, and virus-like particles (VLPs) — require cold-chain storage to maintain potency. Excipients (small-molecule additives) can stabilize these products, but excipient selection is a costly trial-and-error process. Molecular dynamics (MD) simulations can reveal the atomistic mechanisms of excipient–virus interactions, yet even the smallest fully assembled capsids contain **millions of atoms**, requiring massive supercomputing resources for meaningful sampling (see table below).

| Virus | System Size | Performance | Year |
|---|---|---|---|
| STMV | 1.2M atoms | 30 ns/day on 2,016 cores | 2012 |
| HIV-1 | 64M atoms | not reported | 2017 |
| SARS-CoV-2 | 305M atoms | 68 ns/day on 180k cores | 2020 |

**CapSACIN** solves this problem by leveraging icosahedral capsid symmetry to construct minimal surface models focused on a **region of interest (ROI)** embedded in peripheral, context-providing proteins. The resulting models are **20–23% the size of a full capsid** and achieve a **5–18× simulation speed-up** while preserving atomistic resolution and faithfully reproducing both intraprotein and interprotein structure and dynamics.

The method was validated on porcine parvovirus (PPV, PDB: 1K3V) and shown to correctly rank-order experimental excipient stabilization effects (Pearson ρ = −0.974).

---

## The CapSACIN Workflow

The workflow consists of six steps, of which **Steps 1–3 are implemented in this repository**. Steps 4–6 are performed with an external MD engine (GROMACS).

```
┌─────────────────────────────────────────────────────────────┐
│  STEP 1          STEP 2           STEP 3                    │
│  Symmetry-Based  Surface          Position Restraint        │
│  Alignment       Abstraction      Generation                │
│  ─────────────   ─────────────    ──────────────────        │
│  sliceCapsid.py  sliceCapsid.py   genRestraints.py          │
│       │               │                 │                    │
│       ▼               ▼                 ▼                    │
│  Aligned PDB →  Sliced surface  →  GROMACS .itp files       │
│                                                             │
│  ───────────────────────────────────────────────────────    │
│  STEP 4          STEP 5           STEP 6                    │
│  Restrained MD   Unrestrained MD  Nanofragmentation         │
│  (GROMACS)       (GROMACS)        (GROMACS pulling)         │
└─────────────────────────────────────────────────────────────┘
```

### Step 1: Symmetry-Based Alignment

The input capsid (a PDB file) is reoriented so the chosen symmetry axis is perpendicular to the xy-plane (aligned with the z-axis). This is achieved by:

1. Selecting a **reference atom** within the ROI.
2. Finding **two symmetry-related atoms** at the same ROI across the icosahedral capsid (chosen based on fold type: 2-fold, 3-fold, or 5-fold).
3. Computing the **plane normal** defined by the three reference atoms via cross product.
4. Applying **Rodrigues' rotation formula** to align the normal vector with `[0, 0, 1]`.

**Relevant code:** `sliceCapsid.py` lines 55–168, `capsacin/definePlane.py`

### Step 2: Surface Abstraction

The aligned capsid is sliced at a height defined by a **weight parameter ω ∈ [0, 1]**:

- Atoms with `rz < ω · max(z)` are removed.
- ω = 0 retains the full capsid; ω = 1 removes all atoms.
- For PPV, ω = 0.7 provides a 2:1 ratio of peripheral to ROI proteins.

After slicing, a cleanup procedure removes:
- **Fragmented chains** (monomers with fewer residues than expected from the input structure).
- **Broken residues** at the cut boundary (residues with fewer heavy atoms than the canonical count for their amino acid type).

The remaining surface is translated so its minimum z-coordinate is at z = 0, and chain IDs are remapped to a compact character set (A–Z, 0–9, symbols).

**Relevant code:** `sliceCapsid.py` lines 198–310, `capsacin/createDictionary.py`

### Step 3: Position Restraint Generation

To preserve capsid structural integrity during solvent equilibration while allowing the ROI to remain flexible, **scaled position restraints** are applied using an inverse sigmoid function:

```
s(rz,i) = 1 / (1 + exp(rz,i / κ))
K_i = s_i · K_max        (K_max = 1000 kJ mol⁻¹ nm⁻², κ = 0.1)
```

The restraint scheme creates three zones:

| z-range | Restraint Strength | Purpose |
|---|---|---|
| z < 1.5 nm | Full (K = 1000) | Lock peripheral base |
| 1.5–3.0 nm | Scaled via sigmoid | Smooth transition zone |
| z > 3.0 nm | Zero (K = 0) | Fully flexible ROI |

An implicit wall (GROMACS Walls) is placed at z = 0 and z = h_box, preventing excipient diffusion into the capsid interior through periodic boundaries.

**Relevant code:** `genRestraints.py`

### Steps 4 & 5: Restrained and Unrestrained MD Simulations

Performed with **GROMACS** (not included in this repository):

- **Step 4 (Restrained):** 50 ns NVT equilibration at 300 K with capsid position restraints active. Allows solvent to equilibrate while preserving capsid architecture.
- **Step 5 (Unrestrained):** 50–400 ns NVT production with capsid restraints removed (excipient flat-bottom restraints retained). The ROI equilibrates in response to solvent.

Simulation parameters: CHARMM36 force field, TIP3P water, 0.15 M NaCl, PME electrostatics (1.2 nm cutoff), force-switched vdW (1.0–1.2 nm).

### Step 6: Nanofragmentation Simulations

Inspired by Ghaemi et al. (2021), radial pulling forces are applied along vectors from the surface center-of-mass to each monomer's COM. This induces **controlled crack formation** along protein–protein interfaces, enabling assessment of:

- Which symmetry axes are weakest (PPV: 2-fold ≪ 3-fold < 5-fold).
- How excipients modulate interfacial stability (ΔQ_IF).

The fraction of native interfacial contacts, Q_IF, is tracked over the pulling trajectory. The cumulative difference ΔQ_IF between excipient and reference solutions shows **excellent rank-order correlation** with experimental log reduction values (LRV) from thermal stability assays (Pearson ρ = −0.974).

---

## Repository Structure

```
capSACIN/
├── README.md                          # This file
├── LICENSE                            # MIT License
├── desktop/                           # Tauri + React desktop application
│   ├── build_app.sh                   # Complete macOS packaging workflow
│   ├── src/                           # React controls and Mol* viewer
│   └── src-tauri/                     # Rust host, bundle config, and icons
├── env/
│   ├── create-env.sh                  # Conda environment setup script
│   ├── requirements.txt               # Python runtime dependencies
│   └── requirements-dev.txt           # Packaging and test dependencies
└── systemSetup/
    ├── examples.dat                   # Example CLI invocations for various capsids
    ├── sliceCapsid.py                 # ★ Main script: Steps 1 & 2 (alignment + slicing)
    ├── genRestraints.py               # ★ Step 3 (position restraint generation)
    ├── input/                         # 13 bundled capsid PDB structures
    ├── capsacin/                      # Reusable alignment and slicing pipeline
    │   ├── pipeline.py                # Desktop/CLI-compatible computation pipeline
    │   ├── protocol.py                # Sidecar request and result schema
    │   ├── findSymmetryAxes.py        # Global 2/3/5-fold axis search
    │   └── definePlane.py             # Plane and reference-point geometry
    ├── sidecar/                       # JSON-line Python service bundled with the app
    └── tests/                         # Pipeline and sidecar regression tests
```

---

## Module Documentation

### `sliceCapsid.py` — Main Workflow Script

Implements **Steps 1 and 2** of the CapSACIN workflow. This is the primary entry point.

**CLI Arguments:**

| Argument | Type | Default | Description |
|---|---|---|---|
| `--pdb` | str | `9jjh` | PDB filename prefix (reads from `input/{pdb}.pdb`) |
| `--symmetry` | int | `5` | Symmetry axis: `2`, `3`, or `5` (2-fold, 3-fold, or 5-fold) |
| `--refindex` | int[] | `[1058]` | VMD-format atom indices for the reference atom in the ROI |
| `--weight` | float | `0.5` | Slicing fraction ω ∈ [0,1]; 0 = no slice, 0.5 = half, 1 = empty |
| `--plot` | flag | off | Enable 3D matplotlib visualizations for debugging |
| `--auto` | flag | off | Automatically detect the requested symmetry axis and reference atoms |
| `--axis-index` | int | `0` | Candidate axis rank to use; ROI-ranked when `--roi-selection` is set |
| `--roi-selection` | str | `None` | MDAnalysis selection for ROI-aware automatic axis selection |
| `--roi-frame` | int | `0` | MODEL/frame containing the physical ROI copy described by `--roi-selection` |
| `--list-axes` | flag | off | Print candidate axes and exit unless `--auto` is also set |

**Algorithm flow:**
1. Load PDB via `MDAnalysis.Universe`
2. Identify the reference atom (`pointA`) and two symmetry-related atoms (`pointB`, `pointC`)
   - The selection of `pointB` and `pointC` depends on `--symmetry`:
     - 2-fold: nearest and 2nd-nearest same-atom neighbors
     - 3-fold: 2nd- and 3rd-nearest neighbors (k1=1, k2=2)
     - 5-fold: 2nd- and 4th-nearest neighbors (k1=1, k2=3)
3. Compute the plane normal of the three reference atoms (`definePlane.py`)
4. Rotate all coordinates so the normal aligns with `[0, 0, 1]` (Rodrigues' rotation)
5. Center coordinates: COM at origin, min z = 0, min x,y ≥ 0
6. Slice by z-coordinate: keep atoms with `z ≥ ω · max(z)`
7. Clean up broken chains and incomplete residues at the cut boundary
8. Remap chain IDs and save to `output/{pdb}-sliced-sym{symmetry}-w{weight}.pdb`

**ROI-aware automatic axis selection:**

When `--roi-selection` is supplied with `--auto`, CapSACIN ranks all candidate
axes of the requested fold by proximity to the ROI center on `--roi-frame`
(default: frame 0, matching the legacy reference-atom workflow). The selected
axis is oriented toward the ROI before alignment so the retained high-z surface
is centered on the chosen ROI.

Choose an ROI that is local to the expected symmetry feature. If the ROI is
far from every candidate axis, or lies roughly perpendicular to all candidates,
the ranking still returns the closest geometric match, but the result is less
biophysically meaningful and should be inspected with `--list-axes`.

For 2-fold automatic alignment, CapSACIN keeps the legacy convention where
`pointB` comes from the same monomer as `pointA` and is separated by a small
z-perturbation before alignment. This preserves backward compatibility with
the original 2-fold slicing heuristic.

```bash
# Inspect the 5-fold axes closest to residues 250-290 on MODEL/frame 0
python sliceCapsid.py --pdb 1k3v --symmetry 5 --roi-selection "protein and resid 250:290" --list-axes

# Slice using the top ROI-ranked 5-fold axis
python sliceCapsid.py --pdb 1k3v --symmetry 5 --weight 0.7 --auto --roi-selection "protein and resid 250:290"
```

**Output:** A PDB file in `output/` containing the sliced surface model.

---

### `genRestraints.py` — Position Restraint Generator

Implements **Step 3** of the CapSACIN workflow.

**CLI Arguments:**

| Argument | Type | Default | Description |
|---|---|---|---|
| `--pdb` | str | *required* | PDB filename prefix used for `sliceCapsid.py` output |
| `--weight` | float | `1.0` | Global scaling factor for restraint force constants |
| `--symmetry` | int | `None` | Symmetry used for the sliced PDB name; omit only for legacy `output/{pdb}-sliced-w{weight}.pdb` files |

**Algorithm:**
1. Read the sliced PDB from `sliceCapsid.py` output
2. Classify atoms by z-coordinate:
   - **z < 15 Å:** Full restraint (`K = 1000 × weight`)
   - **15 ≤ z < 30 Å:** Scaled via inverted sigmoid: `s = −1/(1+exp(−z_norm/0.1)) + 1`
   - **z ≥ 30 Å:** No restraint (implicitly excluded)
3. Generate per-chain `.itp` files in GROMACS `[ position_restraints ]` format

**Output:** Directory `output/{pdb}-w{weight}-posre/` containing one `.itp` file per chain.

---

### `capsacin/definePlane.py`

Computes the **normal vector** of the plane defined by three 3D points using the cross product `u × v`, where `u` and `v` are vectors from the center-of-mass to two of the points.

**Function:** `definePlane(x, y, z)` → returns the normal vector (length proportional to the area spanned by u,v).

Optionally generates a 3D matplotlib plot showing the points, plane, and normal vector.

---

### `capsacin/formatPDB.py`

Formats a pandas DataFrame's columns to conform to the **fixed-width PDB file format**. Handles:

- Coordinate precision (x, y, z → 8.3f)
- Occupancy and B-factor (→ 6.2f)
- Atom name, residue name, chain ID field widths and justification
- Atom serial number padding

**Function:** `formatPDB(df)` → returns the formatted DataFrame.

---

### `capsacin/createDictionary.py`

Two dictionary utilities:

- **`createDictionary(u)`**: Queries an MDAnalysis Universe for chain A and builds a dictionary mapping each of the 20 standard amino acid three-letter codes to its expected **heavy atom count**. Used to detect broken residues at the slicing boundary.

- **`createChainDictionary(chainValues)`**: Maps a list of old chain IDs to a compact pool of characters (A–Z, 0–9, then symbols: `!@#$%^&-_=+;:'",.<>?/\\|\`~`). The PDB format only reliably supports single-character chain IDs; this remapping ensures output compatibility (and warns if >60 unique chains are present).

---

### `capsacin/alignSymmetry.py` (Legacy)

The original monolithic development script from which `sliceCapsid.py` was refactored. It contains the identical core algorithm but with **hardcoded parameters** instead of CLI arguments. Retained for reference:

- Lines 22–59 contain an extensive commented-out library of `indicesVMD` values for different PDB entries and symmetry types. This serves as a historical reference for atom indices used in various capsid structures.
- Includes additional debugging plots not present in the CLI version.

**Note:** This file is not required for the workflow; `sliceCapsid.py` is the canonical entry point.

---

### `capsacin/mdaCIF.py`

A custom **MDAnalysis reader/parser** for CIF (Crystallographic Information File) format, adapted from [Richard Gowers' implementation](https://github.com/richardjgowers/MDAnalysis_CIFReading). Uses OpenBabel's `pybel` to read crystal structures and unpack symmetry operators (`FillUnitCell`).

- **`CIFReader`**: Reads coordinates and unit cell parameters.
- **`CIFParser`**: Reads atom types, masses, and partial charges.

**Note:** Not currently used by the main workflow. Included for potential future support of crystalline material inputs.

---

## Input PDB Structures

The `input/` directory contains 13 experimentally resolved icosahedral virus capsid structures from the Protein Data Bank:

| PDB ID | Virus / System |
|---|---|
| `1k3v` | **Porcine parvovirus (PPV)** — the primary model system |
| `1dzl` | Additional capsid for validation |
| `1wcd` | Additional capsid for validation |
| `2buk` | Additional capsid for validation |
| `2ztn` | Additional capsid for validation |
| `3r0r` | Additional capsid for validation |
| `3ra2` | Additional capsid for validation |
| `4oq8` | Additional capsid for validation |
| `5cw0` | Additional capsid for validation |
| `6jja` | Additional capsid for validation |
| `8des` | Additional capsid for validation |
| `9clj` | Additional capsid for validation |
| `9jjh` | Additional capsid for validation |

See the paper (Figure S1) for validation of the alignment + slicing step across these diverse capsid architectures.

---

## Command-Line and Development Setup

Desktop users should install the packaged application from the [latest release](https://github.com/Youchengw/capSACIN_modifed/releases/latest). The following setup is for command-line use, testing, or desktop application development.

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/) with Python 3.9
- Node.js 18 or later for desktop frontend development
- Rust and Cargo for native desktop builds

### Python CLI setup

```bash
# 1. Clone the repository
git clone https://github.com/Youchengw/capSACIN_modifed.git
cd capSACIN

# 2. Create and activate the conda environment
cd env
source create-env.sh
conda activate capSACIN

# 3. Verify installation
cd ../systemSetup
python -c "from capsacin import definePlane, formatPDB, createDictionary; print('Ready!')"
```

### Desktop development

```bash
cd desktop
npm ci
npm run check
npm run tauri dev
```

To build the complete Apple Silicon `.app` and `.dmg`, including the Python sidecar and bundled PDB files:

```bash
cd desktop
./build_app.sh
```

The build script requires the `capSACIN` Conda environment, PyInstaller, Node.js, and the Rust toolchain. Generated dependencies, sidecar binaries, application bundles, and Tauri build output are intentionally excluded from Git.

### Dependencies

| Package | Version | Purpose |
|---|---|---|
| `numpy` | 1.26.4 | Array operations, linear algebra |
| `matplotlib` | 3.7.1 | Optional 3D visualization (`--plot` flag) |
| `MDAnalysis` | 2.4.3 | PDB parsing, atom selection, distance computation |
| `tqdm` | 4.65.0 | Progress bars for chain/residue iteration |
| `pandas` | 1.5.3 | Tabular atom data manipulation |
| `pytest` | 8.x | Pipeline and sidecar regression tests |
| `PyInstaller` | 6.x | Standalone Python sidecar packaging |

---

## Command-Line Quick Start

The file `examples.dat` provides ready-to-use command-line examples. Here are the most common use cases:

### 1. Slice a capsid surface model

```bash
cd systemSetup

# PPV 5-fold surface (the primary model from the paper)
python sliceCapsid.py --pdb 1k3v --refindex 959 --symmetry 5 --weight 0.7

# PPV 3-fold surface
python sliceCapsid.py --pdb 1k3v --refindex 3131 --symmetry 3 --weight 0.7

# PPV 2-fold surface
python sliceCapsid.py --pdb 1k3v --refindex 4073 --symmetry 2 --weight 0.7

# Try different slicing weights
python sliceCapsid.py --pdb 1k3v --refindex 959 --symmetry 5 --weight 0.3   # more retained
python sliceCapsid.py --pdb 1k3v --refindex 959 --symmetry 5 --weight 0.55  # moderate

# Other capsids (see examples.dat for more)
python sliceCapsid.py --pdb 9jjh --weight 0.5 --refindex 1058 --symmetry 5
python sliceCapsid.py --pdb 3ra2 --weight 0.5 --refindex 928 --symmetry 5
```

**Output:** `output/{pdb}-sliced-sym{symmetry}-w{weight}.pdb` — the sliced surface model.

### 2. Generate position restraints

```bash
python genRestraints.py --pdb 1k3v --symmetry 5 --weight 0.7
```

**Output:** `output/1k3v-w0.7-posre/` — GROMACS `.itp` files, one per protein chain.

### 3. Running MD simulations (external)

After generating the surface model and restraints, proceed with GROMACS:

```bash
# Step 4: Restrained equilibration (50 ns NVT)
gmx grompp -f restrained.mdp -c surface.pdb -p topol.top -o restrained.tpr
gmx mdrun -deffnm restrained -v

# Step 5: Unrestrained production (50–400 ns NVT)
gmx grompp -f unrestrained.mdp -c restrained.gro -p topol.top -o prod.tpr
gmx mdrun -deffnm prod -v

# Step 6: Nanofragmentation (pulling)
gmx grompp -f pull.mdp -c prod.gro -p topol.top -o pull.tpr
gmx mdrun -deffnm pull -v
```

---

## Key Results from the Paper

Using the CapSACIN workflow with PPV (PDB: 1K3V), the authors demonstrated:

- **Surface models are 20–23% the size** of a fully assembled capsid and achieve **5–18× faster MD sampling** (Table 2 in paper).
- **Intraprotein dynamics** (dynamic cross-correlation) of ROI proteins match the assembled capsid with cos(θ) = 0.964, compared to 0.689 for subunit models lacking peripheral context.
- **Interprotein structure** (sphericity Ψ, monomer COM distances) is preserved within the ROI.
- **Nanofragmentation** identified the 2-fold axis as the weakest interface, followed by the 3-fold, with the 5-fold most resistant — consistent with hierarchical capsid assembly pathways.
- **Excipient rank-ordering** from simulations (ΔQ_IF) matches experimental thermal stability data with **Pearson ρ = −0.974**: sorbitol > trehalose > glutamate > arginine > glycine (most to least stabilizing).

---

## Limitations

- Currently supports only **icosahedral**, nonenveloped viruses with known PDB structures.
- Does not capture **long-range concerted motions** between distant ROIs.
- Nanofragmentation results are **sensitive to pulling parameters** and initial configurations.
- The surface model is **larger** than alternative approaches (e.g., rotational symmetry boundary conditions) but is **MD-engine agnostic**.

---

## Citation

If you use this repository in your research, please cite:

> Jonathan W. P. Zajac, Idris Tohidian, Praveen Muralikrishnan, Caryn L. Heldt, Sarah L. Perry, Sapna Sarupria. "Cracking the Capsid Code: A Computationally Feasible Approach for Investigating Virus–Excipient Interactions in Biologics Design." *J. Chem. Theory Comput.* **2026**, 22, 2635–2651. doi: [10.1021/acs.jctc.5c01810](https://doi.org/10.1021/acs.jctc.5c01810)

```bibtex
@article{Zajac2026,
  title = {Cracking the Capsid Code: A Computationally Feasible Approach for Investigating Virus−Excipient Interactions in Biologics Design},
  author = {Zajac, Jonathan W. P. and Tohidian, Idris and Muralikrishnan, Praveen and Perry, Sarah L. and Heldt, Caryn L. and Sarupria, Sapna},
  journal = {J. Chem. Theory Comput.},
  year = {2026},
  volume = {22},
  pages = {2635--2651},
  doi = {10.1021/acs.jctc.5c01810}
}
```

---

## License

- **Source code and software:** [MIT License](LICENSE)
- **Written and graphical materials:** [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/)

---

## Contact

For questions or feedback, please contact Jonathan Zajac at [zajac028@umn.edu](mailto:zajac028@umn.edu).

The CapSACIN workflow is developed by the [SAMPEL Group](https://sarupria.chem.umn.edu/) (Sarupria Lab) at the University of Minnesota, Twin Cities, in collaboration with the Heldt Lab (Michigan Technological University) and the Perry Lab (University of Massachusetts Amherst).
