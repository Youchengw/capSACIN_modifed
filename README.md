# capSACIN Studio

A desktop interface for capsid slicing, with automatic 2/3/5-fold previews and an integrated Mol* viewer.

Built on [SAMPEL-Group/capSACIN](https://github.com/SAMPEL-Group/capSACIN). For the capSACIN method, scientific background, and citation, see the original repository.

## Download

**v0.1.2 · Apple Silicon · macOS 13+**

- [Download DMG](https://github.com/Youchengw/capSACIN_modifed/releases/download/v0.1.2/capSACIN-Studio-v0.1.2-macOS-arm64.dmg) — open and drag the app into Applications.
- [ZIP, checksums, and release details](https://github.com/Youchengw/capSACIN_modifed/releases/tag/v0.1.2).

Python and 13 example PDBs are bundled. No separate Conda installation is needed. The app is ad-hoc signed and not notarized; if macOS blocks opening it, allow it in **System Settings → Privacy & Security**.

## Use

1. Open a built-in or local PDB. The app automatically prepares all three folds.
2. Switch folds to view cached previews, then adjust the slicing weight.
3. Click **Run capSACIN** to compute the sliced structure.
4. Inspect **Original**, **Sliced**, or **Overlay**, then choose **Save PDB As…**.

Axis or ROI changes require refreshed previews. Manual Reference mode prepares the selected fold. The cache lasts for the currently loaded structure in the current session; startup and rendering still take time. See [release notes](RELEASE_NOTES.md).

Inputs must be multi-MODEL PDBs containing spatial capsid copies with compatible atom ordering, rather than MD trajectory snapshots. Automatic detection does not expand a single asymmetric unit or accept a single-MODEL assembly as equivalent input.

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
