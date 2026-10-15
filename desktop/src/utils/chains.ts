import type { RunSliceResult } from "../sidecar/types.ts";

export function appliedChainIndices(slice: RunSliceResult | null): number[] {
  return (slice?.chain_candidates ?? [])
    .filter((chain) => chain.selected)
    .map((chain) => chain.source_index)
    .sort((a, b) => a - b);
}

export function weightChainIndices(slice: RunSliceResult | null): number[] {
  return (slice?.chain_candidates ?? [])
    .filter((chain) => chain.selected_by_weight)
    .map((chain) => chain.source_index)
    .sort((a, b) => a - b);
}

export function chainSelectionDirty(slice: RunSliceResult | null, selected: number[]): boolean {
  const applied = appliedChainIndices(slice);
  const current = [...new Set(selected)].sort((a, b) => a - b);
  return applied.length !== current.length || applied.some((value, index) => value !== current[index]);
}
