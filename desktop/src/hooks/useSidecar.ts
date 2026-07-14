/**
 * useSidecar hook — manages the Python sidecar lifecycle from React.
 *
 * Wraps the Tauri invoke calls and event listeners for:
 * - inspect_structure
 * - prepare_preview
 * - run_slice
 * - cancellation
 */
import { useCallback, useEffect, useRef } from "react";
import type { UnlistenFn } from "@tauri-apps/api/event";
import type {
  InspectResult,
  PreparePreviewResult,
  RunSliceResult,
  SidecarProgress,
  SidecarResult,
  SliceParams,
} from "../sidecar/types";

interface UseSidecarCallbacks {
  onProgress?: (evt: SidecarProgress) => void;
  onResult?: (evt: SidecarResult) => void;
  onError?: (error: string) => void;
}

export function useSidecar(callbacks?: UseSidecarCallbacks) {
  const unlistenRef = useRef<{ progress?: UnlistenFn; result?: UnlistenFn }>({});

  // Set up event listeners on mount
  useEffect(() => {
    let cancelled = false;

    async function setup() {
      try {
        const { onProgress, onResult } = await import("../sidecar/client");

        const pUnlisten = await onProgress((evt) => {
          if (!cancelled) callbacks?.onProgress?.(evt);
        });

        const rUnlisten = await onResult((evt) => {
          if (!cancelled) {
            if (evt.status === "error" && callbacks?.onError) {
              callbacks.onError(evt.error?.message ?? "Unknown error");
            } else {
              callbacks?.onResult?.(evt);
            }
          }
        });

        unlistenRef.current = { progress: pUnlisten, result: rUnlisten };
      } catch (e) {
        console.warn("[useSidecar] Tauri event listeners not available in browser:", e);
      }
    }

    setup();

    return () => {
      cancelled = true;
      unlistenRef.current.progress?.();
      unlistenRef.current.result?.();
    };
  }, []);

  const inspect = useCallback(async (inputPath: string): Promise<InspectResult | null> => {
    try {
      const { inspectStructure } = await import("../sidecar/client");
      return await inspectStructure(inputPath);
    } catch (e: any) {
      callbacks?.onError?.(e?.message ?? String(e));
      return null;
    }
  }, []);

  const preparePreview = useCallback(async (params: SliceParams): Promise<PreparePreviewResult | null> => {
    try {
      const { preparePreview } = await import("../sidecar/client");
      return await preparePreview(params);
    } catch (e: any) {
      callbacks?.onError?.(e?.message ?? String(e));
      return null;
    }
  }, []);

  const runSlice = useCallback(async (params: SliceParams): Promise<RunSliceResult | null> => {
    try {
      const { runSlice } = await import("../sidecar/client");
      return await runSlice(params);
    } catch (e: any) {
      callbacks?.onError?.(e?.message ?? String(e));
      return null;
    }
  }, []);

  const cancel = useCallback(async () => {
    try {
      const { cancelOperation } = await import("../sidecar/client");
      await cancelOperation();
    } catch (e) {
      console.error("Cancel failed:", e);
    }
  }, []);

  return { inspect, preparePreview, runSlice, cancel };
}
