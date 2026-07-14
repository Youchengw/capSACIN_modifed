import { useAppStore } from "../../store/useAppStore";
import { resolveViewerRoi } from "../../utils/roi";

export function ResultsSummary() {
  const previewResult = useAppStore((s) => s.previewResult);
  const sliceResult = useAppStore((s) => s.sliceResult);
  const inspectResult = useAppStore((s) => s.inspectResult);
  const roiChain = useAppStore((s) => s.roiChain);
  const roiStartResid = useAppStore((s) => s.roiStartResid);
  const roiEndResid = useAppStore((s) => s.roiEndResid);

  if (!previewResult && !sliceResult) return null;

  const stats = sliceResult?.statistics;
  const resolvedRoi = resolveViewerRoi(
    inspectResult,
    roiChain,
    roiStartResid,
    roiEndResid,
  );
  const roiFrames = previewResult?.diagnostics.roi_symmetry_frames;

  return (
    <div className="panel-section">
      <h3>Results</h3>

      {previewResult && !sliceResult && (
        <div style={{ fontSize: "12px", color: "var(--text-secondary)" }}>
          <div>✓ {previewResult.axis_candidates.length} axis candidates found</div>
          <div>✓ Selected axis [{previewResult.selected_axis.map((v) => v.toFixed(3)).join(", ")}]</div>
          {resolvedRoi.chain && Array.isArray(roiFrames) && (
            <div>
              ✓ ROI MODEL {roiFrames.join("/")} · Chain {resolvedRoi.chain} · Residues {resolvedRoi.startResid}–{resolvedRoi.endResid}
            </div>
          )}
          <div style={{ marginTop: "4px" }}>Preview ready — adjust weight, then click Run</div>
        </div>
      )}

      {stats && (
        <div>
          <div className="stats-grid">
            <div className="stat-item">
              <div className="stat-label">Input Atoms</div>
              <div className="stat-value">{stats.input_atoms.toLocaleString()}</div>
            </div>
            <div className="stat-item">
              <div className="stat-label">Output Atoms</div>
              <div className="stat-value">{stats.output_atoms.toLocaleString()}</div>
            </div>
            <div className="stat-item">
              <div className="stat-label">Retained</div>
              <div className="stat-value" style={{ color: "var(--accent)" }}>
                {stats.retained_percentage.toFixed(1)}%
              </div>
            </div>
            <div className="stat-item">
              <div className="stat-label">Output Chains</div>
              <div className="stat-value">{stats.output_chains}</div>
            </div>
            <div className="stat-item">
              <div className="stat-label">Input Residues</div>
              <div className="stat-value">{stats.input_residues.toLocaleString()}</div>
            </div>
            <div className="stat-item">
              <div className="stat-label">Output Residues</div>
              <div className="stat-value">{stats.output_residues.toLocaleString()}</div>
            </div>
          </div>

          {sliceResult && (
            <div style={{ marginTop: "8px", fontSize: "11px", color: "var(--text-secondary)" }}>
              <div>Fold: {useAppStore.getState().symmetry}-fold</div>
              <div>Weight ω: {useAppStore.getState().weight}</div>
              <div>Axis: [{sliceResult.axis.map((v) => v.toFixed(3)).join(", ")}]</div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
