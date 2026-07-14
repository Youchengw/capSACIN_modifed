/** Shared TypeScript types — mirror Python protocol.py dataclasses. */

export type SymmetryFold = 2 | 3 | 5;
export type AxisMode = "auto" | "manual";
export type DisplayMode = "original" | "sliced" | "overlay";
export type OperationType = "inspect_structure" | "prepare_preview" | "run_slice";
export type SidecarStatus = "idle" | "running" | "error" | "cancelled";

// ---- Structure inspection ----

export interface ChainInfo {
  id: string;
  residue_count: number;
  min_resid: number;
  max_resid: number;
  atom_count: number;
}

export interface InspectResult {
  chains: ChainInfo[];
  n_models: number;
  n_chains_per_frame: number;
  total_atoms: number;
  total_residues: number;
  valid_residue_ranges: Record<string, [number, number]>;
  input_path: string;
}

// ---- Axis candidates ----

export interface AxisCandidate {
  axis_index: number;
  axis: [number, number, number];
  score: number;
  ref_index?: number;
  ref_frame?: number;
  chain_id?: string;
  roi_score?: number;
  roi_angle_deg?: number;
}

// ---- Prepare preview ----

export interface PreparePreviewResult {
  axis_candidates: AxisCandidate[];
  selected_axis_index: number;
  selected_axis: [number, number, number];
  aligned_cif_path: string;
  workspace_path: string;
  diagnostics: Record<string, unknown>;
  warnings: string[];
}

// ---- Run slice ----

export interface SliceStatistics {
  input_atoms: number;
  output_atoms: number;
  input_residues: number;
  output_residues: number;
  input_chains: number;
  output_chains: number;
  retained_percentage: number;
}

export interface RunSliceResult {
  output_pdb_path: string;
  sliced_cif_path: string;
  workspace_path: string;
  statistics: SliceStatistics;
  axis: [number, number, number];
  diagnostics: Record<string, unknown>;
  warnings: string[];
}

// ---- Sidecar protocol ----

export interface SidecarProgress {
  type: "progress";
  id: string;
  stage: string;
  fraction: number;
  message: string;
}

export interface SidecarResult {
  type: "result";
  id: string;
  status: "ok" | "error" | "cancelled";
  data?: any;
  error?: { code: string; message: string };
}

export type SidecarEvent = SidecarProgress | SidecarResult;

// ---- Slice parameters ----

export interface SliceParams {
  input_path: string;
  symmetry: SymmetryFold;
  auto: boolean;
  axis_index: number;
  roi_selection: string | null;
  roi_frame: number;
  ref_indices: number[] | null;
  legacy_plane_heuristic: boolean;
  weight: number;
}
