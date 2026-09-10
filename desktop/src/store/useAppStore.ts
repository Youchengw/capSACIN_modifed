import { create } from "zustand";
import type {
  SymmetryFold, AxisMode, DisplayMode, SidecarStatus,
  InspectResult, PreparePreviewResult, RunSliceResult,
  AxisCandidate, SliceStatistics,
  PrepareAllPreviewsResult,
} from "../sidecar/types";

// ---- Built-in PDB list ----

export const BUILTIN_PDBS = [
  "1k3v", "3ra2", "8des", "2ztn", "1dzl", "1wcd",
  "6jja", "9clj", "4oq8", "3r0r", "2buk", "9jjh", "5cw0",
];

// ---- Store shape ----

export interface AppState {
  // Structure
  inputPath: string;
  isBuiltin: boolean;
  structureName: string;
  inputGeneration: number;

  // Parameters
  symmetry: SymmetryFold;
  weight: number;
  axisMode: AxisMode;
  axisIndex: number;
  axisIndices: Record<SymmetryFold, number>;
  roiSelection: string;
  roiEnabled: boolean;
  roiFrame: number;
  roiChain: string;
  roiStartResid: number;
  roiEndResid: number;
  refIndices: string;
  legacyPlane: boolean;

  // Advanced panel
  showAdvanced: boolean;

  // State
  sidecarStatus: SidecarStatus;
  error: string | null;
  stageName: string;
  stageFraction: number;

  // Results
  inspectResult: InspectResult | null;
  previewResult: PreparePreviewResult | null;
  foldPreviews: Partial<Record<SymmetryFold, PreparePreviewResult>>;
  foldErrors: Partial<Record<SymmetryFold, string>>;
  sliceResult: RunSliceResult | null;

  // Viewer
  displayMode: DisplayMode;
  showAxis: boolean;
  showPlane: boolean;
  showRoi: boolean;

  // Actions
  setInputPath: (path: string, isBuiltin: boolean, name: string) => void;
  setSymmetry: (s: SymmetryFold) => void;
  setWeight: (w: number) => void;
  setAxisMode: (m: AxisMode) => void;
  setAxisIndex: (i: number) => void;
  setRoiSelection: (s: string) => void;
  setRoiEnabled: (enabled: boolean) => void;
  setRoiFrame: (frame: number) => void;
  setRoiChain: (chain: string) => void;
  setRoiStartResid: (resid: number) => void;
  setRoiEndResid: (resid: number) => void;
  setRefIndices: (s: string) => void;
  setLegacyPlane: (b: boolean) => void;
  toggleAdvanced: () => void;
  setSidecarStatus: (s: SidecarStatus) => void;
  setProgress: (stage: string, fraction: number) => void;
  setError: (e: string | null) => void;
  setInspectResult: (r: InspectResult) => void;
  setPreviewResult: (r: PreparePreviewResult) => void;
  setAllPreviews: (r: PrepareAllPreviewsResult) => void;
  setSliceResult: (r: RunSliceResult) => void;
  setDisplayMode: (m: DisplayMode) => void;
  toggleAxis: () => void;
  togglePlane: () => void;
  toggleRoi: () => void;
  reset: () => void;
}

const initialState = {
  inputPath: "",
  isBuiltin: false,
  structureName: "",
  inputGeneration: 0,
  symmetry: 5 as SymmetryFold,
  weight: 0.5,
  axisMode: "auto" as AxisMode,
  axisIndex: 0,
  axisIndices: { 2: 0, 3: 0, 5: 0 },
  roiSelection: "",
  roiEnabled: false,
  roiFrame: 0,
  roiChain: "",
  roiStartResid: 0,
  roiEndResid: 0,
  refIndices: "",
  legacyPlane: false,
  showAdvanced: false,
  sidecarStatus: "idle" as SidecarStatus,
  error: null as string | null,
  stageName: "",
  stageFraction: 0,
  inspectResult: null as InspectResult | null,
  previewResult: null as PreparePreviewResult | null,
  foldPreviews: {} as AppState["foldPreviews"],
  foldErrors: {} as AppState["foldErrors"],
  sliceResult: null as RunSliceResult | null,
  displayMode: "original" as DisplayMode,
  showAxis: true,
  showPlane: true,
  showRoi: false,
};

const clearedPreviews = {
  previewResult: null, sliceResult: null, foldPreviews: {}, foldErrors: {}, error: null,
};
const resetAxes = { axisIndex: 0, axisIndices: { 2: 0, 3: 0, 5: 0 } };

export const useAppStore = create<AppState>((set) => ({
  ...initialState,

  setInputPath: (path, isBuiltin, name) =>
    set((state) => ({
      ...clearedPreviews,
      ...resetAxes,
      inputPath: path,
      inputGeneration: state.inputGeneration + 1,
      sidecarStatus: "idle",
      stageName: "",
      stageFraction: 0,
      displayMode: "original",
      isBuiltin,
      structureName: name,
      axisMode: "auto",
      roiEnabled: false,
      roiSelection: "",
      roiFrame: 0,
      roiChain: "",
      roiStartResid: 0,
      roiEndResid: 0,
      refIndices: "",
      inspectResult: null,
      previewResult: null,
      sliceResult: null,
    })),

  setSymmetry: (symmetry) => set((state) => symmetry === state.symmetry ? {} : ({
    symmetry,
    axisIndex: state.axisIndices[symmetry],
    previewResult: state.foldPreviews[symmetry] ?? null,
    sliceResult: null,
    error: state.foldErrors[symmetry] ?? null,
  })),
  setWeight: (weight) => set({ weight, sliceResult: null }),
  setAxisMode: (axisMode) => set({ axisMode, ...clearedPreviews }),
  setAxisIndex: (axisIndex) => set((state) => ({
    axisIndex,
    axisIndices: { ...state.axisIndices, [state.symmetry]: axisIndex },
    previewResult: null, sliceResult: null, error: null,
    foldPreviews: { ...state.foldPreviews, [state.symmetry]: undefined },
    foldErrors: { ...state.foldErrors, [state.symmetry]: undefined },
  })),
  setRoiSelection: (roiSelection) => set((s) => ({ roiSelection, ...(s.roiEnabled ? { ...resetAxes, ...clearedPreviews } : {}) })),
  setRoiEnabled: (roiEnabled) => set({ roiEnabled, ...resetAxes, ...clearedPreviews }),
  setRoiFrame: (roiFrame) => set((s) => ({ roiFrame, ...(s.roiEnabled ? { ...resetAxes, ...clearedPreviews } : {}) })),
  setRoiChain: (roiChain) => set((s) => ({ roiChain, ...(s.roiEnabled ? { ...resetAxes, ...clearedPreviews } : {}) })),
  setRoiStartResid: (roiStartResid) => set((s) => ({ roiStartResid, ...(s.roiEnabled ? { ...resetAxes, ...clearedPreviews } : {}) })),
  setRoiEndResid: (roiEndResid) => set((s) => ({ roiEndResid, ...(s.roiEnabled ? { ...resetAxes, ...clearedPreviews } : {}) })),
  setRefIndices: (refIndices) => set((s) => ({ refIndices, ...(s.axisMode === "manual" ? clearedPreviews : {}) })),
  setLegacyPlane: (legacyPlane) => set({ legacyPlane, ...clearedPreviews }),
  toggleAdvanced: () => set((s) => ({ showAdvanced: !s.showAdvanced })),

  setSidecarStatus: (sidecarStatus) => set({ sidecarStatus }),
  setProgress: (stageName, stageFraction) => set({ stageName, stageFraction }),
  setError: (error) =>
    set((state) => ({
      error,
      sidecarStatus: error ? "error" : state.sidecarStatus,
    })),
  setInspectResult: (inspectResult) => set((state) => {
    const selectedChain = inspectResult.chains.find(
      (chain) => chain.id === state.roiChain,
    ) ?? inspectResult.chains[0];
    return {
      inspectResult,
      roiFrame: Math.min(
        Math.max(0, state.roiFrame),
        Math.max(0, inspectResult.n_models - 1),
      ),
      roiChain: selectedChain?.id ?? "",
      roiStartResid: selectedChain?.min_resid ?? 0,
      roiEndResid: selectedChain?.max_resid ?? 0,
    };
  }),
  setPreviewResult: (previewResult) => set((state) => ({
    previewResult, sliceResult: null, displayMode: "original",
    foldPreviews: { ...state.foldPreviews, [state.symmetry]: previewResult },
    foldErrors: { ...state.foldErrors, [state.symmetry]: undefined },
  })),
  setAllPreviews: (result) => set((state) => ({
    foldPreviews: result.previews,
    foldErrors: result.errors,
    previewResult: result.previews[state.symmetry] ?? null,
    sliceResult: null,
    error: result.errors[state.symmetry] ?? null,
    displayMode: "original",
  })),
  setSliceResult: (sliceResult) => set({ sliceResult }),
  setDisplayMode: (displayMode) => set({ displayMode }),
  toggleAxis: () => set((s) => ({ showAxis: !s.showAxis })),
  togglePlane: () => set((s) => ({ showPlane: !s.showPlane })),
  toggleRoi: () => set((s) => ({ showRoi: !s.showRoi })),
  reset: () => set(initialState),
}));
