import { useMemo, useState } from "react";
import { useAppStore } from "../../store/useAppStore";
import { applyChainSelection } from "../../sidecar/applyChainSelection";
import { chainSelectionDirty, weightChainIndices } from "../../utils/chains";

type SortKey = "z_mean" | "z_min" | "id";

export function ChainSelectionPanel() {
  const slice = useAppStore((state) => state.sliceResult);
  const selected = useAppStore((state) => state.selectedChainIndices);
  const setChainSelected = useAppStore((state) => state.setChainSelected);
  const setSelected = useAppStore((state) => state.setSelectedChainIndices);
  const busy = useAppStore((state) => state.sidecarStatus === "running");
  const [sortKey, setSortKey] = useState<SortKey>("z_mean");
  const candidates = slice?.chain_candidates ?? [];
  const ordered = useMemo(() => [...candidates].sort((a, b) => {
    if (sortKey === "id") return a.id.localeCompare(b.id, undefined, { numeric: true });
    return b[sortKey] - a[sortKey] || a.id.localeCompare(b.id, undefined, { numeric: true });
  }), [candidates, sortKey]);
  if (!slice || !candidates.length) return null;

  const selectedSet = new Set(selected);
  const dirty = chainSelectionDirty(slice, selected);
  return (
    <section className="panel-section chain-selection-panel">
      <h3>Chain selection</h3>
      <p className="chain-note">
        Complete aligned chains are listed by z position. Check chains to add or remove them from the final surface model.
      </p>
      <div className="chain-selection-summary">
        <strong>{selected.length}</strong> / {candidates.length} chains selected
        {dirty && <span>Unapplied changes</span>}
      </div>
      <div className="chain-selection-controls">
        <label>Sort
          <select value={sortKey} disabled={busy} onChange={(event) => setSortKey(event.target.value as SortKey)}>
            <option value="z_mean">Mean z ↓</option>
            <option value="z_min">Minimum z ↓</option>
            <option value="id">Chain ID</option>
          </select>
        </label>
        <div>
          <button className="chain-link" disabled={busy} onClick={() => setSelected(candidates.map((chain) => chain.source_index))}>All</button>
          <button className="chain-link" disabled={busy} onClick={() => setSelected([])}>None</button>
          <button className="chain-link" disabled={busy} onClick={() => setSelected(weightChainIndices(slice))}>Weight result</button>
        </div>
      </div>
      <div className="chain-list">
        {ordered.map((chain) => (
          <label className={`chain-row ${selectedSet.has(chain.source_index) ? "selected" : ""}`} key={chain.source_index}>
            <input type="checkbox" disabled={busy} checked={selectedSet.has(chain.source_index)}
              onChange={(event) => setChainSelected(chain.source_index, event.target.checked)} />
            <span className="chain-id">{chain.id}</span>
            <span className="chain-z">mean {chain.z_mean.toFixed(1)} Å<small>{chain.z_min.toFixed(1)}–{chain.z_max.toFixed(1)} Å</small></span>
            {chain.selected_by_weight && <span className="weight-badge">ω</span>}
          </label>
        ))}
      </div>
      {!selected.length && <p className="chain-note chain-error">Select at least one chain.</p>}
      <button className="btn btn-primary" disabled={busy || !dirty || !selected.length}
        onClick={() => applyChainSelection()}>
        Apply chain selection
      </button>
      <p className="chain-note">ω marks chains retained by the original slicing weight. Applying regenerates the PDB and preview; restraints will use the edited model.</p>
    </section>
  );
}
