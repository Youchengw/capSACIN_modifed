import { useAppStore } from "../store/useAppStore.ts";
import { buildSliceParams } from "../utils/sliceParams.ts";
import type * as Client from "./client";

type PreviewClient = Pick<typeof Client, "prepareAllPreviews" | "preparePreview" | "onProgress">;

/** Used both when opening a PDB and when refreshing edited axis/ROI settings. */
export async function preparePreviews(clientOverride?: PreviewClient) {
  const state = useAppStore.getState();
  if (!state.inputPath || state.sidecarStatus === "running") return;
  const generation = state.inputGeneration;
  const current = () => useAppStore.getState().inputGeneration === generation;
  const params = buildSliceParams(state);
  let unlisten: (() => void) | undefined;
  state.setSidecarStatus("running");
  state.setProgress("Loading PDB and preparing previews…", 0);
  state.setError(null);
  try {
    const client = clientOverride ?? await import("./client");
    unlisten = await client.onProgress((event) => {
      if (current()) state.setProgress(event.stage, event.fraction);
    });
    if (params.auto) {
      const result = await client.prepareAllPreviews({ ...params, axis_indices: state.axisIndices });
      if (!current()) return;
      // Inspection initializes the ROI form only for a newly opened structure.
      if (!state.inspectResult) state.setInspectResult(result.inspection);
      state.setAllPreviews(result);
    } else {
      const result = await client.preparePreview(params);
      if (!current()) return;
      state.setPreviewResult(result);
    }
    state.setSidecarStatus("idle");
  } catch (error) {
    if (!current()) return;
    const message = error instanceof Error ? error.message : String(error);
    if (message.toLowerCase().includes("cancel")) {
      state.setSidecarStatus("cancelled");
      state.setProgress("Cancelled", 0);
    } else {
      state.setError(message);
    }
  } finally {
    unlisten?.();
  }
}
