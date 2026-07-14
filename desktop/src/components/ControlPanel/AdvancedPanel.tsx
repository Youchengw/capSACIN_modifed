import { useAppStore } from "../../store/useAppStore";

export function AdvancedPanel() {
  const showAdvanced = useAppStore((s) => s.showAdvanced);
  const toggleAdvanced = useAppStore((s) => s.toggleAdvanced);
  const axisMode = useAppStore((s) => s.axisMode);
  const setAxisMode = useAppStore((s) => s.setAxisMode);
  const axisIndex = useAppStore((s) => s.axisIndex);
  const setAxisIndex = useAppStore((s) => s.setAxisIndex);
  const roiSelection = useAppStore((s) => s.roiSelection);
  const setRoiSelection = useAppStore((s) => s.setRoiSelection);
  const roiEnabled = useAppStore((s) => s.roiEnabled);
  const setRoiEnabled = useAppStore((s) => s.setRoiEnabled);
  const roiFrame = useAppStore((s) => s.roiFrame);
  const setRoiFrame = useAppStore((s) => s.setRoiFrame);
  const roiChain = useAppStore((s) => s.roiChain);
  const setRoiChain = useAppStore((s) => s.setRoiChain);
  const roiStartResid = useAppStore((s) => s.roiStartResid);
  const setRoiStartResid = useAppStore((s) => s.setRoiStartResid);
  const roiEndResid = useAppStore((s) => s.roiEndResid);
  const setRoiEndResid = useAppStore((s) => s.setRoiEndResid);
  const refIndices = useAppStore((s) => s.refIndices);
  const setRefIndices = useAppStore((s) => s.setRefIndices);
  const legacyPlane = useAppStore((s) => s.legacyPlane);
  const setLegacyPlane = useAppStore((s) => s.setLegacyPlane);
  const previewResult = useAppStore((s) => s.previewResult);
  const inspectResult = useAppStore((s) => s.inspectResult);
  const busy = useAppStore((s) => s.sidecarStatus) === "running";

  const activeRange = inspectResult?.valid_residue_ranges[roiChain];

  return (
    <div className="panel-section">
      <div className="advanced-toggle" onClick={toggleAdvanced}>
        <span style={{ transform: showAdvanced ? "rotate(90deg)" : "none", transition: "transform 0.2s" }}>
          ▸
        </span>
        Advanced Settings
      </div>

      {showAdvanced && (
        <div style={{ marginTop: "8px", display: "flex", flexDirection: "column", gap: "10px" }}>
          {/* Axis Mode */}
          <div>
            <label style={{ fontSize: "11px", color: "var(--text-secondary)" }}>Axis Mode</label>
            <select
              value={axisMode}
              onChange={(e) => setAxisMode(e.target.value as "auto" | "manual")}
              disabled={busy}
            >
              <option value="auto">Automatic</option>
              <option value="manual">Manual Reference</option>
            </select>
          </div>

          {/* Axis Index */}
          {axisMode === "auto" && (
            <div>
              <label style={{ fontSize: "11px", color: "var(--text-secondary)" }}>
                Axis Candidate (rank)
              </label>
              <select
                value={axisIndex}
                onChange={(e) => setAxisIndex(Number(e.target.value))}
                disabled={busy}
              >
                {previewResult
                  ? previewResult.axis_candidates.map((c, i) => (
                      <option key={i} value={i}>
                        Rank {i} (axis {c.axis_index}, score={c.score.toFixed(4)})
                      </option>
                    ))
                  : [0, 1, 2, 3, 4, 5].map((i) => (
                      <option key={i} value={i}>Axis {i}</option>
                    ))}
              </select>
            </div>
          )}

          <div style={{ display: "flex", alignItems: "center", gap: "8px" }}>
            <input
              type="checkbox"
              id="roi-axis-ranking"
              checked={roiEnabled}
              onChange={(e) => setRoiEnabled(e.target.checked)}
              disabled={busy}
            />
            <label
              htmlFor="roi-axis-ranking"
              style={{ fontSize: "11px", color: "var(--text-secondary)", cursor: "pointer" }}
            >
              Use ROI to rank symmetry axes
            </label>
          </div>

          {/* ROI form */}
          <div style={{ display: "flex", gap: "8px" }}>
            <div style={{ flex: 1 }}>
              <label style={{ fontSize: "11px", color: "var(--text-secondary)" }}>ROI Frame</label>
              <input
                type="number"
                min="0"
                max={Math.max(0, (inspectResult?.n_models ?? 1) - 1)}
                value={roiFrame}
                onChange={(e) => setRoiFrame(Number(e.target.value))}
                disabled={busy || !roiEnabled}
              />
            </div>
            <div style={{ flex: 1 }}>
              <label style={{ fontSize: "11px", color: "var(--text-secondary)" }}>Chain</label>
              <select
                value={roiChain}
                onChange={(e) => {
                  const chain = e.target.value;
                  setRoiChain(chain);
                  const range = inspectResult?.valid_residue_ranges[chain];
                  if (range) {
                    setRoiStartResid(range[0]);
                    setRoiEndResid(range[1]);
                  }
                }}
                disabled={busy || !roiEnabled || !inspectResult}
              >
                {(inspectResult?.chains ?? []).map((chain) => (
                  <option key={chain.id} value={chain.id}>{chain.id}</option>
                ))}
              </select>
            </div>
          </div>

          <div style={{ display: "flex", gap: "8px" }}>
            <div style={{ flex: 1 }}>
              <label style={{ fontSize: "11px", color: "var(--text-secondary)" }}>Start residue</label>
              <input
                type="number"
                value={roiStartResid}
                min={activeRange?.[0]}
                max={activeRange?.[1]}
                onChange={(e) => setRoiStartResid(Number(e.target.value))}
                disabled={busy || !roiEnabled || !inspectResult}
              />
            </div>
            <div style={{ flex: 1 }}>
              <label style={{ fontSize: "11px", color: "var(--text-secondary)" }}>End residue</label>
              <input
                type="number"
                value={roiEndResid}
                min={activeRange?.[0]}
                max={activeRange?.[1]}
                onChange={(e) => setRoiEndResid(Number(e.target.value))}
                disabled={busy || !roiEnabled || !inspectResult}
              />
            </div>
          </div>
          {activeRange && (
            <div style={{ fontSize: "10px", color: "var(--text-secondary)" }}>
              Valid range for chain {roiChain}: {activeRange[0]}–{activeRange[1]}
            </div>
          )}

          <div>
            <label style={{ fontSize: "11px", color: "var(--text-secondary)" }}>
              Raw ROI selection (optional override)
            </label>
            <input
              type="text"
              placeholder='e.g. protein and name CA and resid 250:290'
              value={roiSelection}
              onChange={(e) => setRoiSelection(e.target.value)}
              disabled={busy || !roiEnabled}
              style={{ fontSize: "12px" }}
            />
          </div>

          {/* Manual Ref Indices */}
          {axisMode === "manual" && (
            <div>
              <label style={{ fontSize: "11px", color: "var(--text-secondary)" }}>
                Reference Atom Indices
              </label>
              <input
                type="text"
                placeholder="e.g. 959, 1058, 1200"
                value={refIndices}
                onChange={(e) => setRefIndices(e.target.value)}
                disabled={busy}
              />
            </div>
          )}

          {/* Legacy Plane */}
          <div style={{ display: "flex", alignItems: "center", gap: "8px" }}>
            <input
              type="checkbox"
              id="legacy-plane"
              checked={legacyPlane}
              onChange={(e) => setLegacyPlane(e.target.checked)}
              disabled={busy}
            />
            <label htmlFor="legacy-plane" style={{ fontSize: "11px", color: "var(--text-secondary)", cursor: "pointer" }}>
              Legacy plane heuristic (3-point + z-perturbation)
            </label>
          </div>
        </div>
      )}
    </div>
  );
}
