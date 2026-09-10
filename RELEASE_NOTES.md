# capSACIN Studio v0.1.2

Opening a PDB now automatically prepares the 2-fold, 3-fold, and 5-fold previews. Switching between prepared folds reuses their results and updates the 3D view without running Prepare again.

## Changes

- Load the structure once per preparation request and share its global symmetry-axis search and atom data across all three folds.
- Generate candidate axes and scores without computing unused per-candidate reference-atom diagnostics. The CLI retains its full reference diagnostics.
- Reuse prepared previews when switching folds or adjusting weight. Remember each fold's selected axis candidate.
- Show preparation progress and readiness for each fold. A failed fold remains unavailable while successful folds can still be selected.
- Preserve manual reference preparation for the selected fold. Changes to axis selection or ROI ranking settings require refreshed previews.
- Clear obsolete views when changing structures and ignore late preparation responses for an earlier PDB.
- Speed up mmCIF generation while preserving MODEL/chain identities and coordinate formatting.
- Fix the packaged engine entry point so multiprocessing helper processes do not start duplicate preview servers.
- Give only the Python engine the library-loading entitlement needed for its bundled libraries, retaining Hardened Runtime on both the engine and app. Exercise the signed engine during release builds.
- Explicitly select the existing 13 bundled PDBs, keeping additional local structures out of installers.

## Validation

- Python pipeline and sidecar regression tests, including all three folds and ROI-aware alignment.
- Exact comparison of batch-preview coordinates with independently aligned coordinates at the existing output precision.
- Frontend tests for cached fold switching, parameter invalidation, partial failure, cancellation, and late responses.
- TypeScript and native build checks, packaged application startup, and packaged engine execution.

## Performance and limitations

A local test using 1k3v prepared all three folds in about 7.8 seconds after the packaged engine was ready. The previous single 5-fold preparation took about 30.4 seconds in a source-level measurement. These are individual local measurements, not a general speed guarantee, and exclude 3D rendering.

The engine's first startup remains slow in the tested packaged build (about 34 seconds in one run). Prepared-fold switching does not restart the engine or repeat preparation. The viewer still needs time to parse and render a fold's structure.

## Downloads and requirements

- `capSACIN-Studio-v0.1.2-macOS-arm64.dmg` — recommended installer.
- `capSACIN-Studio-v0.1.2-macOS-arm64.zip` — zipped application bundle.
- `SHA256SUMS.txt` — checksums, continuing the existing release format.

Requires an Apple Silicon Mac and macOS 13 or later. Python and the 13 example structures are bundled; a separate Conda environment is not required.

This build is ad-hoc signed and is not notarized with an Apple Developer ID. If macOS blocks the first launch, right-click the app and choose **Open**, or allow it from **System Settings → Privacy & Security**.
