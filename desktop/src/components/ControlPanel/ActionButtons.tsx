import { useAppStore } from "../../store/useAppStore";
import { buildSliceParams } from "../../utils/sliceParams";
import { preparePreviews } from "../../sidecar/preparePreviews";

export function ActionButtons() {
  const inputPath = useAppStore((s) => s.inputPath);
  const axisMode = useAppStore((s) => s.axisMode);
  const sidecarStatus = useAppStore((s) => s.sidecarStatus);
  const previewResult = useAppStore((s) => s.previewResult);
  const sliceResult = useAppStore((s) => s.sliceResult);
  const setSidecarStatus = useAppStore((s) => s.setSidecarStatus);
  const setProgress = useAppStore((s) => s.setProgress);
  const setError = useAppStore((s) => s.setError);
  const setSliceResult = useAppStore((s) => s.setSliceResult);

  const busy = sidecarStatus === "running";
  const hasStructure = !!inputPath;
  const hasPreview = !!previewResult;

  const handleOperationError = (error: unknown) => {
    const message = error instanceof Error ? error.message : String(error);
    if (message.toLowerCase().includes("cancel")) {
      setError(null);
      setSidecarStatus("cancelled");
      setProgress("Cancelled", 0);
      return;
    }
    setError(message);
  };

  const handleRunSlice = async () => {
    if (!hasStructure) return;
    let unlisten: (() => void) | undefined;
    try {
      setSidecarStatus("running");
      setError(null);

      const { runSlice, onProgress } = await import("../../sidecar/client");
      unlisten = await onProgress((evt) => {
        setProgress(evt.stage, evt.fraction);
      });

      const params = buildSliceParams(useAppStore.getState());
      const result = await runSlice(params);
      setSliceResult(result);
      setSidecarStatus("idle");
    } catch (error) {
      handleOperationError(error);
    } finally {
      unlisten?.();
    }
  };

  const handleCancel = async () => {
    try {
      const { cancelOperation } = await import("../../sidecar/client");
      await cancelOperation();
      setSidecarStatus("cancelled");
      setProgress("Cancelled", 0);
    } catch (e) {
      console.error("Cancel failed:", e);
    }
  };

  const handleSave = async () => {
    if (!sliceResult) return;
    try {
      const { savePdbDialog } = await import("../../sidecar/client");
      await savePdbDialog(sliceResult.output_pdb_path);
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    }
  };

  return (
    <div className="panel-section">
      <div style={{ display: "flex", gap: "8px", flexWrap: "wrap" }}>
        {!busy && (
          <button
            className="btn btn-primary"
            disabled={!hasStructure}
            onClick={() => preparePreviews()}
            style={{ flex: 1 }}
          >
            {axisMode === "auto" ? (hasPreview ? "Refresh previews" : "Prepare all folds") : "Prepare Preview"}
          </button>
        )}

        {busy && (
          <button className="btn btn-danger" onClick={handleCancel} style={{ flex: 1 }}>
            Cancel
          </button>
        )}

        {!busy && hasPreview && (
          <button
            className="btn btn-success"
            onClick={handleRunSlice}
            style={{ flex: 1 }}
          >
            Run capSACIN
          </button>
        )}

        {!busy && sliceResult && (
          <button
            className="btn btn-outline"
            onClick={handleSave}
            style={{ flex: 1 }}
          >
            Save PDB As…
          </button>
        )}
      </div>
    </div>
  );
}
