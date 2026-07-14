import { create } from "zustand";
import type {
  SymmetryFold, AxisMode, DisplayMode, SidecarStatus,
  InspectResult, PreparePreviewResult, RunSliceResult,
  AxisCandidate, SliceStatistics,
} from "../sidecar/types";

// ---- Built-in PDB list ----

export const BUILTIN_PDBS = [
  "1k3v", "3ra2", "8des", "2ztn", "1dzl", "1wcd",
  "6jja", "9clj", "4oq8", "3r0r", "2buk", "9jjh", "5cw0",
];

// ---- Store shape ----

interface AppState {
  // Structure
  inputPath: string;
  isBuiltin: boolean;
  structureName: string;

  // Parameters
  symmetry: SymmetryFold;
  weight: number;
  axisMode: AxisMode;
  axisIndex: number;
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
  symmetry: 5 as SymmetryFold,
  weight: 0.5,
  axisMode: "auto" as AxisMode,
  axisIndex: 0,
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
  sliceResult: null as RunSliceResult | null,
  displayMode: "original" as DisplayMode,
  showAxis: true,
  showPlane: true,
  showRoi: false,
};

export const useAppStore = create<AppState>((set) => ({
  ...initialState,

  setInputPath: (path, isBuiltin, name) =>
    set({
      inputPath: path,
      isBuiltin,
      structureName: name,
      axisIndex: 0,
      roiSelection: "",
      roiFrame: 0,
      roiChain: "",
      roiStartResid: 0,
      roiEndResid: 0,
      refIndices: "",
      inspectResult: null,
      previewResult: null,
      sliceResult: null,
    }),

  setSymmetry: (symmetry) => set({ symmetry, axisIndex: 0, previewResult: null, sliceResult: null }),
  setWeight: (weight) => set({ weight, sliceResult: null }),
  setAxisMode: (axisMode) => set({ axisMode, previewResult: null, sliceResult: null }),
  setAxisIndex: (axisIndex) => set({ axisIndex, previewResult: null, sliceResult: null }),
  setRoiSelection: (roiSelection) => set({ roiSelection, axisIndex: 0, previewResult: null, sliceResult: null }),
  setRoiEnabled: (roiEnabled) => set({ roiEnabled, axisIndex: 0, previewResult: null, sliceResult: null }),
  setRoiFrame: (roiFrame) => set({ roiFrame, axisIndex: 0, previewResult: null, sliceResult: null }),
  setRoiChain: (roiChain) => set({ roiChain, axisIndex: 0, previewResult: null, sliceResult: null }),
  setRoiStartResid: (roiStartResid) => set({ roiStartResid, axisIndex: 0, previewResult: null, sliceResult: null }),
  setRoiEndResid: (roiEndResid) => set({ roiEndResid, axisIndex: 0, previewResult: null, sliceResult: null }),
  setRefIndices: (refIndices) => set({ refIndices, previewResult: null, sliceResult: null }),
  setLegacyPlane: (legacyPlane) => set({ legacyPlane, previewResult: null, sliceResult: null }),
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
  setPreviewResult: (previewResult) => set({ previewResult, displayMode: "original" }),
  setSliceResult: (sliceResult) => set({ sliceResult }),
  setDisplayMode: (displayMode) => set({ displayMode }),
  toggleAxis: () => set((s) => ({ showAxis: !s.showAxis })),
  togglePlane: () => set((s) => ({ showPlane: !s.showPlane })),
  toggleRoi: () => set((s) => ({ showRoi: !s.showRoi })),
  reset: () => set(initialState),
}));
