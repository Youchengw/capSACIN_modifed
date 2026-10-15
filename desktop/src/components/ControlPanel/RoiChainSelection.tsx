import { useMemo } from "react";
import { useAppStore } from "../../store/useAppStore";
import { defaultRoiChainIndices, resolveRoiFrames } from "../../utils/roi";

export function RoiChainSelection() {
  const slice = useAppStore((state) => state.sliceResult);
  const preview = useAppStore((state) => state.previewResult);
  const roiFrame = useAppStore((state) => state.roiFrame);
  const roiChain = useAppStore((state) => state.roiChain);
  const selected = useAppStore((state) => state.roiChainIndices);
  const setChainSelected = useAppStore((state) => state.setRoiChainSelected);
  const setSelected = useAppStore((state) => state.setRoiChainIndices);
  const busy = useAppStore((state) => state.sidecarStatus === "running");
  const candidates = useMemo(
    () => (slice?.chain_candidates ?? [])
      .filter((chain) => chain.selected)
      .sort((a, b) => b.z_mean - a.z_mean || a.id.localeCompare(b.id, undefined, { numeric: true })),
    [slice],
  );
  if (!slice || !candidates.length) return null;

  const selectedSet = new Set(selected);
  const defaults = defaultRoiChainIndices(
    slice,
    resolveRoiFrames(preview?.diagnostics, roiFrame),
    roiChain,
  );
  return (
    <div className="roi-chain-selection">
      <h4>ROI chains</h4>
      <p className="restraint-note">
        Select any retained chains to highlight as the ROI. For example, a 2-fold surface can use 2, 4, or 6 ROI chains.
      </p>
      <div className="chain-selection-summary">
        <strong>{selected.length}</strong> / {candidates.length} ROI chains
      </div>
      <div className="roi-chain-actions">
        <button className="chain-link" disabled={busy} onClick={() => setSelected(defaults)}>Symmetry default</button>
        <button className="chain-link" disabled={busy} onClick={() => setSelected(candidates.map((chain) => chain.source_index))}>All</button>
        <button className="chain-link" disabled={busy} onClick={() => setSelected([])}>None</button>
      </div>
      <div className="chain-list roi-chain-list">
        {candidates.map((chain) => (
          <label className={`chain-row ${selectedSet.has(chain.source_index) ? "selected" : ""}`} key={chain.source_index}>
            <input type="checkbox" disabled={busy} checked={selectedSet.has(chain.source_index)}
              onChange={(event) => setChainSelected(chain.source_index, event.target.checked)} />
            <span className="chain-id">{chain.id}</span>
            <span className="chain-z">mean {chain.z_mean.toFixed(1)} Å<small>{chain.z_min.toFixed(1)}–{chain.z_max.toFixed(1)} Å</small></span>
            {defaults.includes(chain.source_index) && <span className="roi-default-badge" title="Default chain for this symmetry axis">ROI</span>}
          </label>
        ))}
      </div>
      <p className="restraint-note">
        ROI chain selection controls the red overlay only. Restraint strengths still follow the three z regions above.
      </p>
    </div>
  );
}
