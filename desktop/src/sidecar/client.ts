/**
 * Sidecar client — communicates with the Python sidecar process
 * through Tauri commands (invoke).
 */

import { invoke } from "@tauri-apps/api/core";
import { listen, type UnlistenFn } from "@tauri-apps/api/event";
import type {
  InspectResult, PreparePreviewResult, RunSliceResult,
  SidecarProgress, SidecarResult, SliceParams,
} from "./types";

// ---- Event listeners ----

export function onProgress(cb: (evt: SidecarProgress) => void): Promise<UnlistenFn> {
  return listen<SidecarProgress>("sidecar:progress", (e) => cb(e.payload));
}

export function onResult(cb: (evt: SidecarResult) => void): Promise<UnlistenFn> {
  return listen<SidecarResult>("sidecar:result", (e) => cb(e.payload));
}

// ---- Operations ----

export async function inspectStructure(inputPath: string): Promise<InspectResult> {
  return invoke<InspectResult>("inspect_structure", { inputPath });
}

export async function preparePreview(params: SliceParams): Promise<PreparePreviewResult> {
  return invoke<PreparePreviewResult>("prepare_preview", { params });
}

export async function runSlice(params: SliceParams): Promise<RunSliceResult> {
  return invoke<RunSliceResult>("run_slice", { params });
}

export async function cancelOperation(): Promise<void> {
  return invoke("cancel_operation");
}

// ---- File dialogs ----

export async function openPdbDialog(): Promise<string | null> {
  return invoke<string | null>("open_pdb_dialog");
}

export async function savePdbDialog(sourcePath: string): Promise<string | null> {
  return invoke<string | null>("save_pdb_dialog", { sourcePath });
}

export async function listBuiltinStructures(): Promise<string[]> {
  return invoke<string[]>("list_builtin_structures");
}

export async function getStructurePath(name: string): Promise<string> {
  return invoke<string>("get_structure_path", { name });
}

export async function readViewerFile(path: string): Promise<string> {
  return invoke<string>("read_text_file", { path });
}
