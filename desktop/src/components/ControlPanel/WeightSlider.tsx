import { useAppStore } from "../../store/useAppStore";

export function WeightSlider() {
  const weight = useAppStore((s) => s.weight);
  const setWeight = useAppStore((s) => s.setWeight);
  const busy = useAppStore((s) => s.sidecarStatus) === "running";

  return (
    <div className="panel-section">
      <h3>Slicing Weight ω</h3>
      <div style={{ display: "flex", alignItems: "center", gap: "12px" }}>
        <input
          type="range"
          className="weight-slider"
          min="0"
          max="1"
          step="0.01"
          value={weight}
          onChange={(e) => setWeight(parseFloat(e.target.value))}
          disabled={busy}
          style={{ flex: 1 }}
        />
        <input
          type="number"
          min="0"
          max="1"
          step="0.01"
          value={weight}
          onChange={(e) => setWeight(Math.min(1, Math.max(0, parseFloat(e.target.value) || 0)))}
          disabled={busy}
          style={{ width: "60px", textAlign: "center" }}
        />
      </div>
      <p style={{ fontSize: "10px", color: "var(--text-secondary)", marginTop: "4px" }}>
        Higher weight removes more structure — drag to preview, click Run to compute
      </p>
    </div>
  );
}
