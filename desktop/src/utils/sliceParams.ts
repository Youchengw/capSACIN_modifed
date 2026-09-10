import type { AppState } from "../store/useAppStore";
import type { SliceParams } from "../sidecar/types";

export function buildSliceParams(state: AppState): SliceParams {
  const formRoi = state.roiStartResid > 0 && state.roiEndResid >= state.roiStartResid
    ? `protein${state.roiChain ? ` and chainid ${state.roiChain}` : ""} and resid ${state.roiStartResid}:${state.roiEndResid}`
    : "";
  return {
    input_path: state.inputPath,
    symmetry: state.symmetry,
    auto: state.axisMode === "auto",
    axis_index: state.axisIndex,
    roi_selection: state.roiEnabled ? (state.roiSelection.trim() || formRoi || null) : null,
    roi_frame: state.roiFrame,
    ref_indices: state.refIndices ? state.refIndices.split(",").map(Number).filter((n) => !isNaN(n)) : null,
    legacy_plane_heuristic: state.legacyPlane,
    weight: state.weight,
  };
}
