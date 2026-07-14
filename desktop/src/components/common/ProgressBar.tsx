interface Props {
  stage: string;
  fraction: number;
  onCancel?: () => void;
}

const STAGE_LABELS: Record<string, string> = {
  load: "Loading PDB…",
  detect_axis: "Detecting symmetry axis…",
  plane_points: "Computing plane points…",
  extract_pdb: "Extracting atom data…",
  align: "Aligning capsid…",
  format: "Formatting coordinates…",
  slice: "Slicing by weight…",
  cleanup: "Cleaning broken residues…",
  write: "Writing output…",
  inspect: "Inspecting structure…",
  prepare: "Preparing preview…",
};

export function ProgressBar({ stage, fraction, onCancel }: Props) {
  const pct = Math.round(fraction * 100);
  const label = STAGE_LABELS[stage] ?? stage;

  return (
    <div>
      <div style={{ display: "flex", justifyContent: "space-between", marginBottom: "4px" }}>
        <span style={{ fontSize: "11px", color: "var(--text-secondary)" }}>{label}</span>
        <span style={{ fontSize: "11px", color: "var(--accent)", fontWeight: 600 }}>{pct}%</span>
      </div>
      <div className="progress-bar">
        <div
          className="progress-bar-fill"
          style={{ width: `${pct}%` }}
        />
      </div>
    </div>
  );
}
