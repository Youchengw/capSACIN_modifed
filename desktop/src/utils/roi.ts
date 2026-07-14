import type { InspectResult } from "../sidecar/types";

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
