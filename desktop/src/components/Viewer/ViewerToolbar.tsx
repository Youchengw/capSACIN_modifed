import { useAppStore } from "../../store/useAppStore";
import type { DisplayMode } from "../../sidecar/types";

const MODES: { value: DisplayMode; label: string }[] = [
  { value: "original", label: "Original" },
  { value: "sliced", label: "Sliced" },
  { value: "overlay", label: "Overlay" },
];

export function ViewerToolbar() {
  const displayMode = useAppStore((s) => s.displayMode);
  const setDisplayMode = useAppStore((s) => s.setDisplayMode);
  const showAxis = useAppStore((s) => s.showAxis);
  const showPlane = useAppStore((s) => s.showPlane);
  const showRoi = useAppStore((s) => s.showRoi);
  const toggleAxis = useAppStore((s) => s.toggleAxis);
  const togglePlane = useAppStore((s) => s.togglePlane);
  const toggleRoi = useAppStore((s) => s.toggleRoi);
  const previewResult = useAppStore((s) => s.previewResult);
  const sliceResult = useAppStore((s) => s.sliceResult);

  const hasData = !!previewResult || !!sliceResult;

  if (!hasData) return null;

  return (
    <div className="viewer-toolbar">
      {/* Display mode */}
      {MODES.map((m) => (
        <button
          key={m.value}
          className={displayMode === m.value ? "active" : ""}
          onClick={() => setDisplayMode(m.value)}
          title={
            m.value !== "original" && !sliceResult
              ? `${m.label} preview based on the current weight; Run capSACIN for the final result`
              : m.label
          }
        >
          {m.label}
        </button>
      ))}

      <div style={{ width: "1px", background: "rgba(255,255,255,0.2)", margin: "0 4px" }} />

      {/* Toggles */}
      <button className={showAxis ? "active" : ""} onClick={toggleAxis}>
        Axis
      </button>
      <button className={showPlane ? "active" : ""} onClick={togglePlane}>
        Plane
      </button>
      <button className={showRoi ? "active" : ""} onClick={toggleRoi}>
        ROI
      </button>

      <div style={{ width: "1px", background: "rgba(255,255,255,0.2)", margin: "0 4px" }} />

      <button onClick={() => window.dispatchEvent(new Event("capsacin:reset-camera"))}>
        Reset
      </button>
      <button onClick={() => window.dispatchEvent(new Event("capsacin:screenshot"))}>
        📷
      </button>
    </div>
  );
}
