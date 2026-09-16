import { useAppStore } from "../store/useAppStore.ts";
import type * as Client from "./client";

type RestraintClient = Pick<typeof Client, "generateRestraints" | "onProgress">;

export async function generateRestraints(clientOverride?: RestraintClient) {
  const state = useAppStore.getState();
  if (!state.sliceResult || state.sidecarStatus === "running") return;
  const current = () => {
    const now = useAppStore.getState();
    return now.sliceResult === state.sliceResult && now.restraintSettings === state.restraintSettings;
  };
  let unlisten: (() => void) | undefined;
  state.setSidecarStatus("running");
  state.setProgress("Generating restraints…", 0);
  state.setError(null);
  state.setRestraintResult(null);
  try {
    const client = clientOverride ?? await import("./client");
    unlisten = await client.onProgress((event) => {
      if (current()) state.setProgress(event.stage, event.fraction);
    });
    const result = await client.generateRestraints({
      input_path: state.sliceResult.output_pdb_path, ...state.restraintSettings,
    });
    if (!current()) return;
    state.setRestraintResult(result);
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
    // A parameter edit can invalidate this result without starting a new job.
    const now = useAppStore.getState();
    if (!current() && now.inputGeneration === state.inputGeneration && now.sidecarStatus === "running") {
      state.setSidecarStatus("idle");
    }
  }
}
