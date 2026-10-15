# capSACIN Studio

A desktop interface for capsid slicing, with automatic 2/3/5-fold previews and an integrated Mol* viewer.

Built on [SAMPEL-Group/capSACIN](https://github.com/SAMPEL-Group/capSACIN). For the capSACIN method, scientific background, and citation, see the original repository.

## Download

**v0.1.3 · Apple Silicon · macOS 13+**

- [Download DMG](https://github.com/Youchengw/capSACIN_modifed/releases/download/v0.1.3/capSACIN-Studio-v0.1.3-macOS-arm64.dmg) — open and drag the app into Applications.
- [ZIP, checksums, and release details](https://github.com/Youchengw/capSACIN_modifed/releases/tag/v0.1.3).

Python and 13 example PDBs are bundled. No separate Conda installation is needed.

The app is ad-hoc signed and has not been notarized by Apple. If macOS blocks it on first launch, follow [Apple's guide to opening an app from an unknown developer](https://support.apple.com/en-gb/guide/mac-help/mh40616/mac).

## Use

1. Open a built-in or local PDB. The app automatically prepares all three folds.
2. Switch folds to view cached previews, then adjust the slicing weight.
3. Click **Run capSACIN** to compute the sliced structure.
4. In **Chain selection**, review complete chains sorted by z, add or remove chains,
   and click **Apply chain selection** when needed.
5. Inspect **Original**, **Sliced**, or **Overlay**, then choose **Save PDB As…**.

Axis or ROI changes require refreshed previews. Manual Reference mode prepares the selected fold. The cache lasts for the currently loaded structure in the current session; startup and rendering still take time. See [release notes](RELEASE_NOTES.md).

Inputs must be multi-MODEL PDBs containing spatial capsid copies with compatible atom ordering, rather than MD trajectory snapshots. Automatic detection does not expand a single asymmetric unit or accept a single-MODEL assembly as equivalent input.

### Chain selection

After slicing, Studio lists every complete aligned chain, including chains outside
the weight-based result. Chain identities use `M<MODEL>_<original chain>` so they
remain traceable to the input. The list is sorted by mean z by default and also
shows each chain's z range; it can instead be sorted by minimum z or chain ID.
The `ω` badge marks chains retained by the slicing weight. Select any combination
and apply it to regenerate the final PDB and sliced preview without repeating axis
detection or alignment. Restraint generation then uses that edited PDB.

### Restraint regions

After slicing, open **Restraint generation** and click **Generate & preview restraints**.
The **ROI chains** list defaults to the local n-fold copies and can be changed to
any retained-chain combination, such as 2, 4, or 6 chains for a 2-fold surface.
This selection controls the red ROI overlay; restraint strengths remain defined
by the three z regions.
The **Restraints** view colors atoms in three regions using the final sliced PDB:

- **Full restraint** (pink): z < 15 Å.
- **Transition** (yellow): 15 ≤ z < 30 Å, with an inverted sigmoid decay.
- **Free** (green): z ≥ 30 Å; these atoms are omitted from the restraint files.

Two translucent planes mark the region boundaries (15 and 30 Å by default).
Use **Boundaries** in the viewer toolbar to show or hide them. The planes follow
the generated restraint settings and can be toggled without resetting the view.
Use **ROI** to overlay the selected region in red, including in the Restraints
view. Toggling the highlight preserves the camera; turn it off to see the
underlying restraint colors. Region counts and exported restraints are unchanged.

The slice bottom is z = 0 Å. Region boundaries and restraint strength are editable;
the default strength scale is 1, giving full K = 1000 kJ mol⁻¹ nm⁻² in all three
directions. Strength is independent of slicing weight. The sigmoid retains the
CLI's normalization over the observed transition atom heights. A transition with
only one distinct height needs adjusted boundaries because that normalization
is undefined; empty regions are allowed.

**Download restraints (.zip)** saves one GROMACS `.itp` per chain plus a README
with settings and indexing details. Atom numbers start at 1 within each chain in
the sliced PDB order; they must match the corresponding molecule topology after
any atom additions, deletions, or reordering. Changing slice or restraint settings
clears the previous export until regeneration. The command-line generator shares
the same calculation and retains its existing filenames and arguments.

For slices above 99,999 atoms, PDB serials wrap to five digits, matching the
MDAnalysis writer convention. Atom order and chain IDs determine restraint
indices; viewer mmCIF atom IDs remain unique. Restraint generation also reads
older capSACIN slices whose serials overflowed into the `ATOM` record field.

## CLI and development

### Python setup

```bash
git clone https://github.com/Youchengw/capSACIN_modifed.git capSACIN
cd capSACIN/env
source create-env.sh
conda activate capSACIN
cd ../systemSetup

# Example; run CLI commands from systemSetup/
python sliceCapsid.py --pdb 1k3v --symmetry 5 --weight 0.7 --auto
python sliceCapsid.py --help
```

The environment uses Python 3.9. Dependency lists are in [env/](env/); additional CLI examples are in [examples.dat](systemSetup/examples.dat).

### Use with a coding agent

The repository includes a portable [capSACIN CLI skill](skills/capsacin-cli/SKILL.md)
for input inspection, axis listing, automatic/ROI-aware slicing, and manual
reference selection. Set up Python as above; the skill runs the existing engine
and does not require the desktop app.

Codex discovers it through the repository's `.agents/skills/capsacin-cli` link
when working in this checkout ([skill discovery documentation](https://learn.chatgpt.com/docs/build-skills#where-codex-loads-local-skills)).
For example, ask:

```text
Use $capsacin-cli to slice the bundled 1k3v PDB automatically around a fivefold
axis with weight 0.7. Save the PDB and report the retained atom and chain counts.
```

With another coding agent, point it to `skills/capsacin-cli/SKILL.md` and ask it
to follow those instructions using its terminal tools. Agents supporting the
Agent Skills format can also install/copy the entire `skills/capsacin-cli`
folder into their skill directory. When using a copied skill outside this
checkout, supply the repository and Python paths via `--repo`/`--python` or
`CAPSACIN_REPO`/`CAPSACIN_PYTHON`.

The [launcher](skills/capsacin-cli/scripts/run_capsacin.py) accepts an external
PDB path and creates a new run directory containing the staged input, `run.log`,
and `output/`. Native CLI parameters are forwarded unchanged; the launcher
handles the input prefix and working directory. See the
[CLI reference](skills/capsacin-cli/references/cli.md) for command examples,
input limitations, and optional restraint generation.

### Desktop development

Requires an Apple Silicon Mac, Rust/Cargo, and Node.js 22.18+ (22.x) or 24+ with npm. Tested with Node.js 26.4.0. After Python setup, build the engine first; from the repository root:

```bash
conda activate capSACIN
cd systemSetup
python -m PyInstaller sidecar/capsacin_sidecar.spec --noconfirm
mkdir -p ../desktop/src-tauri/binaries
install -m 755 dist/capsacin-sidecar ../desktop/src-tauri/binaries/capsacin-sidecar-aarch64-apple-darwin
cd ../desktop
npm ci
npm run tauri dev
```

Rebuild the engine after Python changes. To generate signed app, DMG, and ZIP packages, run `./build_app.sh` from `desktop/`; outputs go to `releases/v<version>/`.

Checks:

```bash
# From systemSetup/, with capSACIN activated
python -m pytest tests -q

# From desktop/
npm run check
npm test
```

Local literature-review notes, `docs/`, and build outputs remain outside version control. The 13 bundled PDBs are explicitly listed in [tauri.conf.json](desktop/src-tauri/tauri.conf.json).

Source code: [MIT](LICENSE). Written and graphical materials: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
