#!/usr/bin/env bash
set -euo pipefail
# capSACIN Studio — Build script for macOS Apple Silicon
#
# Prerequisites:
#   - Rust toolchain (rustc, cargo)
#   - Node.js >= 18 + npm
#   - Python 3.9+ with conda env "capSACIN"
#   - PyInstaller installed in that environment

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
SYSTEM_SETUP="$PROJECT_ROOT/systemSetup"
TARGET_TRIPLE="aarch64-apple-darwin"
SIDECAR_NAME="capsacin-sidecar-$TARGET_TRIPLE"
SIDECAR_DIST="/tmp/capsacin-sidecar-dist"
export PYINSTALLER_CONFIG_DIR="${PYINSTALLER_CONFIG_DIR:-/tmp/capsacin-pyinstaller-config}"

if [[ -n "${CAPSACIN_PYTHON:-}" ]]; then
    PYTHON=("$CAPSACIN_PYTHON")
elif command -v conda >/dev/null 2>&1; then
    PYTHON=(conda run --no-capture-output -n capSACIN python)
elif [[ -x /opt/anaconda3/bin/conda ]]; then
    PYTHON=(/opt/anaconda3/bin/conda run --no-capture-output -n capSACIN python)
else
    PYTHON=(python3)
fi

export PATH="$HOME/.cargo/bin:$PATH"

echo "=== capSACIN Studio Build ==="
echo ""

# Step 1: Build Python sidecar
echo "[1/6] Building Python sidecar with PyInstaller..."
cd "$SYSTEM_SETUP"
"${PYTHON[@]}" -c "import PyInstaller" >/dev/null || {
    echo "PyInstaller is missing. Install it in the capSACIN environment first."
    exit 1
}
mkdir -p "$SCRIPT_DIR/src-tauri/binaries" "$SIDECAR_DIST"
"${PYTHON[@]}" -m PyInstaller sidecar/capsacin_sidecar.spec \
    --clean --noconfirm \
    --distpath "$SIDECAR_DIST" \
    --workpath /tmp/capsacin-pyinstaller
install -m 755 "$SIDECAR_DIST/capsacin-sidecar" \
    "$SCRIPT_DIR/src-tauri/binaries/$SIDECAR_NAME"
echo "  → Sidecar built: $SIDECAR_NAME"

# Step 2: Validate the release's explicitly selected built-in PDBs. Local
# structures can remain in input/ without silently entering the installer.
echo "[2/6] Checking bundled PDB resources..."
PDB_COUNT=$("${PYTHON[@]}" - "$SCRIPT_DIR/src-tauri/tauri.conf.json" <<'PY'
import json
import sys
from pathlib import Path

config_path = Path(sys.argv[1])
resources = json.loads(config_path.read_text())["bundle"]["resources"]
for source in resources:
    path = (config_path.parent / source).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Bundled PDB not found: {path}")
print(len(resources))
PY
)
if [[ "$PDB_COUNT" -ne 13 ]]; then
    echo "Expected 13 built-in PDB files, found $PDB_COUNT."
    exit 1
fi
echo "  → $PDB_COUNT PDBs will be bundled from systemSetup/input"

# Step 3: Install frontend dependencies
echo "[3/6] Installing npm dependencies..."
cd "$SCRIPT_DIR"
npm ci
echo "  → Dependencies installed"

# Step 4: Check frontend. Tauri's beforeBuildCommand builds it in step 5.
echo "[4/6] Checking frontend..."
npm run check
echo "  → Frontend check passed"

# Step 5: Build Tauri app
echo "[5/6] Building Tauri app bundle..."
cd "$SCRIPT_DIR"
npm run tauri build -- --target "$TARGET_TRIPLE" --bundles app --ci

echo "[6/6] Signing, exercising and packaging the application..."
APP="$SCRIPT_DIR/src-tauri/target/$TARGET_TRIPLE/release/bundle/macos/capSACIN Studio.app"
# Tauri applies the same entitlements to every binary. Sign the Python engine
# separately so only it can load PyInstaller's extracted ad-hoc libraries.
codesign --force --sign - --options runtime \
    --entitlements "$SCRIPT_DIR/src-tauri/sidecar.entitlements.plist" \
    "$APP/Contents/MacOS/capsacin-sidecar"
# Do not use --deep when signing: preserve the engine's separate entitlements.
codesign --force --sign - --options runtime "$APP"
codesign --verify --deep --strict "$APP"
"${PYTHON[@]}" "$SCRIPT_DIR/check_packaged_app.py" "$APP"

VERSION=$("${PYTHON[@]}" -c 'import json,sys; print(json.load(open(sys.argv[1]))["version"])' \
    "$SCRIPT_DIR/src-tauri/tauri.conf.json")
RELEASE_DIR="$PROJECT_ROOT/releases/v$VERSION"
ASSET="capSACIN-Studio-v$VERSION-macOS-arm64"
mkdir -p "$RELEASE_DIR"
ditto "$APP" "$RELEASE_DIR/capSACIN Studio.app"
STAGING=$(mktemp -d /tmp/capsacin-dmg.XXXXXX)
trap 'rm -rf "$STAGING"' EXIT
ditto "$APP" "$STAGING/capSACIN Studio.app"
ln -s /Applications "$STAGING/Applications"
# Package the already signed app directly. Re-running Tauri's DMG bundler
# would replace the engine signature and its dedicated entitlements.
hdiutil create -ov -volname "capSACIN Studio" -srcfolder "$STAGING" \
    -format UDZO "$RELEASE_DIR/$ASSET.dmg"
hdiutil verify "$RELEASE_DIR/$ASSET.dmg"
ditto -c -k --sequesterRsrc --keepParent "$APP" "$RELEASE_DIR/$ASSET.zip"
cp "$PROJECT_ROOT/RELEASE_NOTES.md" "$RELEASE_DIR/RELEASE_NOTES.md"
(cd "$RELEASE_DIR" && shasum -a 256 "$ASSET.dmg" "$ASSET.zip" > SHA256SUMS.txt)

echo ""
echo "=== Build complete ==="
echo "App, DMG, ZIP and checksums: $RELEASE_DIR"
