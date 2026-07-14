import { useAppStore } from "../../store/useAppStore";
import type { SliceParams } from "../../sidecar/types";

export function ActionButtons() {
  const inputPath = useAppStore((s) => s.inputPath);
  const symmetry = useAppStore((s) => s.symmetry);
  const weight = useAppStore((s) => s.weight);
  const axisMode = useAppStore((s) => s.axisMode);
  const axisIndex = useAppStore((s) => s.axisIndex);
  const roiSelection = useAppStore((s) => s.roiSelection);
  const roiEnabled = useAppStore((s) => s.roiEnabled);
  const roiFrame = useAppStore((s) => s.roiFrame);
  const roiChain = useAppStore((s) => s.roiChain);
  const roiStartResid = useAppStore((s) => s.roiStartResid);
  const roiEndResid = useAppStore((s) => s.roiEndResid);
  const refIndices = useAppStore((s) => s.refIndices);
  const legacyPlane = useAppStore((s) => s.legacyPlane);
  const sidecarStatus = useAppStore((s) => s.sidecarStatus);
  const previewResult = useAppStore((s) => s.previewResult);
  const sliceResult = useAppStore((s) => s.sliceResult);
  const setSidecarStatus = useAppStore((s) => s.setSidecarStatus);
  const setProgress = useAppStore((s) => s.setProgress);
  const setError = useAppStore((s) => s.setError);
  const setPreviewResult = useAppStore((s) => s.setPreviewResult);
  const setSliceResult = useAppStore((s) => s.setSliceResult);

  const busy = sidecarStatus === "running";
  const hasStructure = !!inputPath;
  const hasPreview = !!previewResult;

  const formRoiSelection = roiStartResid > 0 && roiEndResid >= roiStartResid
    ? `protein${roiChain ? ` and chainid ${roiChain}` : ""} and resid ${roiStartResid}:${roiEndResid}`
    : "";

  const buildParams = (): SliceParams => ({
    input_path: inputPath,
    symmetry,
    auto: axisMode === "auto",
    axis_index: axisIndex,
    roi_selection: roiEnabled ? (roiSelection.trim() || formRoiSelection || null) : null,
    roi_frame: roiFrame,
    ref_indices: refIndices ? refIndices.split(",").map(Number).filter((n) => !isNaN(n)) : null,
    legacy_plane_heuristic: legacyPlane,
    weight,
  });

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

  const handlePreparePreview = async () => {
    if (!hasStructure) return;
    let unlisten: (() => void) | undefined;
    try {
      setSidecarStatus("running");
      setError(null);

      const { preparePreview, onProgress } = await import("../../sidecar/client");
      unlisten = await onProgress((evt) => {
        setProgress(evt.stage, evt.fraction);
      });

      const params = buildParams();
      const result = await preparePreview(params);
      setPreviewResult(result);
      setSidecarStatus("idle");
    } catch (error) {
      handleOperationError(error);
    } finally {
      unlisten?.();
    }
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

      const params = buildParams();
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
            onClick={handlePreparePreview}
            style={{ flex: 1 }}
          >
            Prepare Preview
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
