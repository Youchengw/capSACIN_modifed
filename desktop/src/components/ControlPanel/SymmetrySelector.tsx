import type { SymmetryFold } from "../../sidecar/types";
import { useAppStore } from "../../store/useAppStore";

const FOLDS: { value: SymmetryFold; label: string; color: string }[] = [
  { value: 2, label: "2-fold", color: "#22c55e" },
  { value: 3, label: "3-fold", color: "#3b82f6" },
  { value: 5, label: "5-fold", color: "#eab308" },
];

export function SymmetrySelector() {
  const symmetry = useAppStore((s) => s.symmetry);
  const setSymmetry = useAppStore((s) => s.setSymmetry);
  const busy = useAppStore((s) => s.sidecarStatus) === "running";
  const previews = useAppStore((s) => s.foldPreviews);
  const errors = useAppStore((s) => s.foldErrors);

  return (
    <div className="panel-section">
      <h3>Symmetry Fold</h3>
      <div className="fold-radio">
        {FOLDS.map((f) => (
          <button
            key={f.value}
            className={`fold-btn ${symmetry === f.value ? "selected" : ""}`}
            disabled={busy}
            title={errors[f.value] ?? (previews[f.value] ? "Preview ready" : "Preview not prepared")}
            onClick={() => setSymmetry(f.value)}
            style={symmetry === f.value ? { background: f.color, borderColor: f.color } : {}}
          >
            {f.label}{previews[f.value] ? " ✓" : errors[f.value] ? " !" : ""}
          </button>
        ))}
      </div>
    </div>
  );
}
