import { useAppStore } from "../store/useAppStore.ts";
import type * as Client from "./client";

type ChainSelectionClient = Pick<typeof Client, "applyChainSelection" | "onProgress">;

export async function applyChainSelection(clientOverride?: ChainSelectionClient) {
  const state = useAppStore.getState();
  if (!state.sliceResult || state.sidecarStatus === "running" || !state.selectedChainIndices.length) return;
  const source = state.sliceResult.chain_source_path;
  const selected = [...state.selectedChainIndices].sort((a, b) => a - b);
  const current = () => {
    const now = useAppStore.getState();
    return now.sliceResult?.chain_source_path === source
      && now.selectedChainIndices.length === selected.length
      && [...now.selectedChainIndices].sort((a, b) => a - b)
        .every((value, index) => value === selected[index]);
  };
  let unlisten: (() => void) | undefined;
  state.setSidecarStatus("running");
  state.setProgress("Applying chain selection…", 0);
  state.setError(null);
  state.setRestraintResult(null);
  try {
    const client = clientOverride ?? await import("./client");
    unlisten = await client.onProgress((event: any) => {
      if (current()) state.setProgress(event.stage, event.fraction);
    });
    const result = await client.applyChainSelection({
      chain_source_path: source,
      selected_chain_indices: selected,
    });
    if (!current()) return;
    state.setSliceResult(result);
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
    const now = useAppStore.getState();
    if (!current() && now.sidecarStatus === "running") state.setSidecarStatus("idle");
  }
}
