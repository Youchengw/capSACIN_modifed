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

# Step 2: Copy bundled PDBs
echo "[2/6] Checking bundled PDB resources..."
PDB_COUNT=$(find "$SYSTEM_SETUP/input" -maxdepth 1 -name '*.pdb' -type f | wc -l | tr -d ' ')
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

# Step 4: Build frontend
echo "[4/6] Building frontend..."
npm run build
echo "  → Frontend built"

# Step 5: Build Tauri app
echo "[5/6] Building Tauri app bundle..."
cd "$SCRIPT_DIR"
npm run tauri build -- --target "$TARGET_TRIPLE"

echo "[6/6] Running packaged-app startup smoke test..."
APP_EXEC="$SCRIPT_DIR/src-tauri/target/$TARGET_TRIPLE/release/bundle/macos/capSACIN Studio.app/Contents/MacOS/capsacin-studio"
CAPSACIN_STARTUP_CHECK=1 "$APP_EXEC"
echo "  → Tauri and plugins initialized successfully"

echo ""
echo "=== Build complete ==="
echo "App:  src-tauri/target/$TARGET_TRIPLE/release/bundle/macos/"
echo "DMG:  src-tauri/target/$TARGET_TRIPLE/release/bundle/dmg/"
