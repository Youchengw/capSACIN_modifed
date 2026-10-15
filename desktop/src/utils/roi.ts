import type { InspectResult, RunSliceResult } from "../sidecar/types";

export interface ResolvedViewerRoi {
  chain: string;
  startResid: number;
  endResid: number;
}

/** Resolve the ROI against the structure that is actually loaded. */
export function resolveViewerRoi(
  inspectResult: InspectResult | null,
  requestedChain: string,
  requestedStart: number,
  requestedEnd: number,
): ResolvedViewerRoi {
  if (!inspectResult || inspectResult.chains.length === 0) {
    return {
      chain: requestedChain,
      startResid: requestedStart,
      endResid: requestedEnd,
    };
  }

  const requested = inspectResult.chains.find(
    (chain) => chain.id === requestedChain,
  );
  const selected = requested ?? inspectResult.chains[0];
  const validRange = inspectResult.valid_residue_ranges[selected.id]
    ?? [selected.min_resid, selected.max_resid];
  const requestedRangeIsValid = Boolean(requested)
    && requestedStart >= validRange[0]
    && requestedEnd <= validRange[1]
    && requestedEnd >= requestedStart;

  return {
    chain: selected.id,
    startResid: requestedRangeIsValid ? requestedStart : validRange[0],
    endResid: requestedRangeIsValid ? requestedEnd : validRange[1],
  };
}

export function resolveRoiFrames(
  diagnostics: Record<string, unknown> | undefined,
  fallbackFrame: number,
): number[] {
  const frames = diagnostics?.roi_symmetry_frames;
  if (!Array.isArray(frames)) return [fallbackFrame];
  const valid = frames.filter(
    (frame): frame is number => Number.isInteger(frame) && frame >= 0,
  );
  return valid.length > 0 ? valid : [fallbackFrame];
}

export function defaultRoiChainIndices(
  slice: RunSliceResult | null,
  frames: number[],
  originalChainId: string,
): number[] {
  if (!slice) return [];
  const frameSet = new Set(frames);
  const applied = (slice.chain_candidates ?? []).filter(
    (chain) => chain.selected && frameSet.has(chain.frame),
  );
  const exact = applied.filter(
    (chain) => chain.original_chain_id === originalChainId,
  );
  return (exact.length > 0 ? exact : applied)
    .map((chain) => chain.source_index)
    .sort((a, b) => a - b);
}

export function selectedRoiChainIds(
  slice: RunSliceResult | null,
  selectedIndices: number[],
): string[] {
  if (!slice) return [];
  const selected = new Set(selectedIndices);
  return (slice.chain_candidates ?? [])
    .filter((chain) => chain.selected && selected.has(chain.source_index))
    .sort((a, b) => b.z_mean - a.z_mean || a.id.localeCompare(b.id, undefined, { numeric: true }))
    .map((chain) => chain.id);
}
