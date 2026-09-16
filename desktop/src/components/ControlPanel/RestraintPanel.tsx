import { useState } from "react";
import { useAppStore } from "../../store/useAppStore";
import { generateRestraints } from "../../sidecar/generateRestraints";
import { RESTRAINT_REGIONS } from "../../utils/restraints";

export function RestraintPanel() {
  const slice = useAppStore((s) => s.sliceResult);
  const result = useAppStore((s) => s.restraintResult);
  const settings = useAppStore((s) => s.restraintSettings);
  const setSettings = useAppStore((s) => s.setRestraintSettings);
  const busy = useAppStore((s) => s.sidecarStatus) === "running";
  const [saved, setSaved] = useState<{ archive: string; path: string } | null>(null);
  if (!slice) return null;
  const valid = Object.values(settings).every(Number.isFinite)
    && settings.full_until >= 0 && settings.free_from > settings.full_until && settings.strength >= 0;
  const ranges = {
    full: `z < ${settings.full_until} Å`,
    transition: `${settings.full_until} ≤ z < ${settings.free_from} Å`,
    free: `z ≥ ${settings.free_from} Å`,
  };
  const save = async () => {
    if (!result) return;
    try {
      const { saveRestraintsDialog } = await import("../../sidecar/client");
      const path = await saveRestraintsDialog(result.archive_path);
      if (path) setSaved({ archive: result.archive_path, path });
    } catch (error) {
      useAppStore.getState().setError(error instanceof Error ? error.message : String(error));
    }
  };
  return (
    <section className="panel-section restraint-panel">
      <h3>Restraint generation</h3>
      <p className="restraint-note">Define three regions along the sliced structure’s z axis. Its bottom is at z = 0 Å.</p>
      <div className="restraint-settings">
        <label>Full until (Å)
          <input aria-label="Full restraint boundary" type="number" min="0" step="1" disabled={busy}
            value={Number.isFinite(settings.full_until) ? settings.full_until : ""} onChange={(e) => setSettings({ full_until: e.target.valueAsNumber })} />
        </label>
        <label>Free from (Å)
          <input aria-label="Free region boundary" type="number" min="0" step="1" disabled={busy}
            value={Number.isFinite(settings.free_from) ? settings.free_from : ""} onChange={(e) => setSettings({ free_from: e.target.valueAsNumber })} />
        </label>
      </div>
      <label>Restraint strength scale
        <input aria-label="Restraint strength scale" type="number" min="0" step="0.1" disabled={busy}
          value={Number.isFinite(settings.strength) ? settings.strength : ""} onChange={(e) => setSettings({ strength: e.target.valueAsNumber })} />
      </label>
      <p className="restraint-note">Full K = {Number.isFinite(settings.strength) ? (1000 * settings.strength).toLocaleString() : "—"} kJ mol⁻¹ nm⁻². Independent of slice weight ω.</p>
      <div className="restraint-regions">
        {RESTRAINT_REGIONS.map((region) => (
          <div className="restraint-region" key={region.key}>
            <span className="region-dot" style={{ background: region.color }} />
            <div><strong>{region.label}</strong><small>{ranges[region.key]}</small></div>
            {result && <span>{result.region_counts[region.key].toLocaleString()} atoms</span>}
          </div>
        ))}
      </div>
      {!valid && <p className="restraint-note" role="alert">Use 0 ≤ full boundary &lt; free boundary and a nonnegative strength.</p>}
      {slice.statistics.output_atoms === 0 && <p className="restraint-note">The slice is empty. Adjust the slice settings first.</p>}
      <button className="btn btn-primary" disabled={busy || !valid || slice.statistics.output_atoms === 0}
        onClick={() => generateRestraints()}>
        {result ? "Regenerate restraints" : "Generate & preview restraints"}
      </button>
      {result && <>
        <p className="restraint-note">{result.files.length} chain files · z = {result.z_range[0]}–{result.z_range[1]} Å</p>
        <button className="btn btn-outline" disabled={busy} onClick={save}>Download restraints (.zip)</button>
        <p className="restraint-note">Includes one GROMACS .itp per chain and settings. Atom indices follow each chain’s order in the sliced PDB.</p>
        {saved?.archive === result.archive_path && <p className="restraint-note" role="status">Saved to {saved.path}</p>}
      </>}
    </section>
  );
}
