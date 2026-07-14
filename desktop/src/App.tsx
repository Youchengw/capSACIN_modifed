import { useAppStore } from "./store/useAppStore";
import { StructureSelector } from "./components/ControlPanel/StructureSelector";
import { SymmetrySelector } from "./components/ControlPanel/SymmetrySelector";
import { WeightSlider } from "./components/ControlPanel/WeightSlider";
import { ActionButtons } from "./components/ControlPanel/ActionButtons";
import { ResultsSummary } from "./components/ControlPanel/ResultsSummary";
import { AdvancedPanel } from "./components/ControlPanel/AdvancedPanel";
import { ProgressBar } from "./components/common/ProgressBar";
import { ErrorBanner } from "./components/common/ErrorBanner";
import { MolStarViewer } from "./components/Viewer/MolStarViewer";
import { ViewerToolbar } from "./components/Viewer/ViewerToolbar";

export default function App() {
  const error = useAppStore((s) => s.error);
  const sidecarStatus = useAppStore((s) => s.sidecarStatus);
  const stageName = useAppStore((s) => s.stageName);
  const stageFraction = useAppStore((s) => s.stageFraction);
  const setSidecarStatus = useAppStore((s) => s.setSidecarStatus);
  const setError = useAppStore((s) => s.setError);

  const cancelCurrentOperation = async () => {
    try {
      const { cancelOperation } = await import("./sidecar/client");
      await cancelOperation();
      setSidecarStatus("cancelled");
    } catch (cancelError: any) {
      setError(cancelError?.message ?? String(cancelError));
    }
  };

  return (
    <div className="app-layout">
      {/* Left Control Panel */}
      <aside className="control-panel">
        <div className="panel-section" style={{ padding: "16px" }}>
          <h1 style={{ fontSize: "18px", fontWeight: 700, color: "var(--text-primary)" }}>
            capSACIN Studio
          </h1>
          <p style={{ fontSize: "11px", color: "var(--text-secondary)", marginTop: "2px" }}>
            Capsid Surface Abstraction &amp; Nanofragmentation
          </p>
        </div>

        <StructureSelector />
        <SymmetrySelector />
        <WeightSlider />

        <ActionButtons />

        {sidecarStatus === "running" && (
          <div className="panel-section">
            <ProgressBar
              stage={stageName}
              fraction={stageFraction}
              onCancel={cancelCurrentOperation}
            />
          </div>
        )}

        {error && (
          <div className="panel-section">
            <ErrorBanner message={error} />
          </div>
        )}

        <ResultsSummary />
        <AdvancedPanel />
      </aside>

      {/* Right 3D Viewer */}
      <main className="viewer-area">
        <MolStarViewer />
        <ViewerToolbar />
      </main>
    </div>
  );
}
