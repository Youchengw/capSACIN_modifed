import { useEffect, useRef, useState } from "react";
import { createPluginUI } from "molstar/lib/mol-plugin-ui";
import { renderReact18 } from "molstar/lib/mol-plugin-ui/react18";
import { DefaultPluginUISpec } from "molstar/lib/mol-plugin-ui/spec";
import { MolScriptBuilder as MS } from "molstar/lib/mol-script/language/builder";
import { Color } from "molstar/lib/mol-util/color";
import type { PluginUIContext } from "molstar/lib/mol-plugin-ui/context";

import { useAppStore } from "../../store/useAppStore";
import { readViewerFile } from "../../sidecar/client";
import { resolveViewerRoi } from "../../utils/roi";

type Bounds = { min: [number, number, number]; max: [number, number, number] };

const COLORS = {
  retained: Color(0x55d6be),
  removed: Color(0x526071),
  sliced: Color(0x59a8ff),
  roi: Color(0xff5d8f),
  axis: Color(0xffcf56),
  plane: Color(0x6ee7f2),
};

export function MolStarViewer() {
  const containerRef = useRef<HTMLDivElement>(null);
  const pluginRef = useRef<PluginUIContext | null>(null);
  const cacheRef = useRef(new Map<string, string>());
  const rebuildSerial = useRef(0);
  const rebuildQueueRef = useRef<Promise<void>>(Promise.resolve());
  const [viewerReady, setViewerReady] = useState(false);
  const [viewerMessage, setViewerMessage] = useState("Initializing 3D viewer…");

  const previewResult = useAppStore((state) => state.previewResult);
  const inputGeneration = useAppStore((state) => state.inputGeneration);
  const inputPath = useAppStore((state) => state.inputPath);
  const preparing = useAppStore((state) => state.sidecarStatus) === "running";
  const sliceResult = useAppStore((state) => state.sliceResult);
  const displayMode = useAppStore((state) => state.displayMode);
  const weight = useAppStore((state) => state.weight);
  const showAxis = useAppStore((state) => state.showAxis);
  const showPlane = useAppStore((state) => state.showPlane);
  const showRoi = useAppStore((state) => state.showRoi);
  const roiFrame = useAppStore((state) => state.roiFrame);
  const roiChain = useAppStore((state) => state.roiChain);
  const roiStartResid = useAppStore((state) => state.roiStartResid);
  const roiEndResid = useAppStore((state) => state.roiEndResid);
  const inspectResult = useAppStore((state) => state.inspectResult);
  const setError = useAppStore((state) => state.setError);
  const resolvedRoi = resolveViewerRoi(
    inspectResult,
    roiChain,
    roiStartResid,
    roiEndResid,
  );

  useEffect(() => { cacheRef.current.clear(); }, [inputGeneration]);

  useEffect(() => {
    let disposed = false;
    let plugin: PluginUIContext | null = null;

    async function initialize() {
      if (!containerRef.current) return;
      plugin = await createPluginUI({
        target: containerRef.current,
        render: renderReact18,
        spec: {
          ...DefaultPluginUISpec(),
          layout: {
            initial: { isExpanded: false, showControls: false },
          },
          components: {
            controls: { left: "none", right: "none", top: "none", bottom: "none" },
            remoteState: "none",
          },
        },
      });
      if (disposed) {
        plugin.dispose();
        return;
      }
      pluginRef.current = plugin;
      setViewerReady(true);
      setViewerMessage("");
    }

    initialize().catch((error) => {
      const message = `Could not initialize Mol*: ${error instanceof Error ? error.message : String(error)}`;
      setViewerMessage(message);
      setError(message);
    });

    const resetCamera = () => pluginRef.current?.canvas3d?.requestCameraReset();
    const saveScreenshot = () => {
      pluginRef.current?.helpers.viewportScreenshot?.download("capsacin-view.png");
    };
    window.addEventListener("capsacin:reset-camera", resetCamera);
    window.addEventListener("capsacin:screenshot", saveScreenshot);

    return () => {
      disposed = true;
      window.removeEventListener("capsacin:reset-camera", resetCamera);
      window.removeEventListener("capsacin:screenshot", saveScreenshot);
      pluginRef.current = null;
      plugin?.dispose();
    };
  }, [setError]);

  useEffect(() => {
    if (!viewerReady || !pluginRef.current) return;
    const serial = ++rebuildSerial.current;
    const timeout = window.setTimeout(() => {
      // Mol* state-tree mutations are asynchronous and are not safe to run in
      // parallel. Serialize rebuilds so an older Preview cannot add its ROI
      // after a newer rebuild has already cleared the plugin.
      rebuildQueueRef.current = rebuildQueueRef.current
        .catch(() => undefined)
        .then(async () => {
          const plugin = pluginRef.current;
          if (!plugin || serial !== rebuildSerial.current) return;
          setViewerMessage("Rendering capsid…");
          try {
            await plugin.clear();
            if (!previewResult && !sliceResult) {
              setViewerMessage("");
              return;
            }
            let originalStructure: any = null;
            let slicedStructure: any = null;
            const needsOriginal =
              displayMode === "original"
              || displayMode === "overlay"
              || (displayMode === "sliced" && !sliceResult);

            if (needsOriginal) {
              if (!previewResult?.aligned_cif_path) {
                throw new Error("Aligned preview is not available.");
              }
              const text = await cachedViewerFile(previewResult.aligned_cif_path, cacheRef.current);
              originalStructure = await loadMmcif(plugin, text, "Aligned capsid");
              if (displayMode === "original") {
                await addWeightPreview(plugin, originalStructure, weight);
              } else if (displayMode === "overlay") {
                await addAllRepresentation(plugin, originalStructure, "overlay-original", COLORS.removed, 0.18);
              }
            }

            if (displayMode === "sliced" || displayMode === "overlay") {
              if (sliceResult?.sliced_cif_path) {
                const text = await cachedViewerFile(sliceResult.sliced_cif_path, cacheRef.current);
                slicedStructure = await loadMmcif(plugin, text, "Sliced capsid");
                await addAllRepresentation(plugin, slicedStructure, "sliced", COLORS.sliced, 0.95);
              } else {
                if (!originalStructure) {
                  throw new Error("Aligned preview is not available.");
                }
                await addRetainedRepresentation(
                  plugin,
                  originalStructure,
                  weight,
                  "estimated-slice",
                  COLORS.sliced,
                );
              }
            }

            const roiStructure = displayMode === "original"
              ? originalStructure
              : (sliceResult ? slicedStructure : originalStructure);
            if (
              showRoi
              && roiStructure
              && resolvedRoi.chain
              && resolvedRoi.startResid > 0
              && resolvedRoi.endResid >= resolvedRoi.startResid
            ) {
              const roiFrames = getRoiFrames(previewResult?.diagnostics, roiFrame);
              await addRoiRepresentation(
                plugin,
                roiStructure,
                roiFrames.map((frame) => `M${frame}_${resolvedRoi.chain}`),
                resolvedRoi.startResid,
                resolvedRoi.endResid,
                displayMode !== "original" && !sliceResult ? weight : undefined,
              );
            }

            const bounds = getBounds(previewResult?.diagnostics);
            if (bounds && (showAxis || showPlane)) {
              const markerCif = buildMarkerCif(bounds, weight, showAxis, showPlane);
              const markerStructure = await loadMmcif(plugin, markerCif, "Guides");
              if (showAxis) await addAsymRepresentation(plugin, markerStructure, "AXIS", "axis", COLORS.axis, 0.95, 1.5);
              if (showPlane) await addAsymRepresentation(plugin, markerStructure, "PLANE", "plane", COLORS.plane, 0.28, 1.0);
            }

            if (serial === rebuildSerial.current) {
              setViewerMessage("");
              plugin.canvas3d?.requestCameraReset();
            }
          } catch (error) {
            if (serial !== rebuildSerial.current) return;
            const message = `3D rendering failed: ${error instanceof Error ? error.message : String(error)}`;
            setViewerMessage(message);
            setError(message);
          }
        });
    }, 180);

    return () => window.clearTimeout(timeout);
  }, [
    viewerReady,
    inputGeneration,
    previewResult,
    sliceResult,
    displayMode,
    weight,
    showAxis,
    showPlane,
    showRoi,
    roiFrame,
    roiChain,
    roiStartResid,
    roiEndResid,
    inspectResult,
    setError,
  ]);

  const hasStructure = Boolean(previewResult || sliceResult);

  return (
    <div className="molstar-shell">
      <div ref={containerRef} className="molstar-host" />
      {!hasStructure && viewerReady && (
        <div className="viewer-empty-state">
          <div className="viewer-empty-icon">◌</div>
          <div>{!inputPath ? "Select a PDB to begin" : preparing ? "Preparing symmetry previews…" : "No preview for the current settings"}</div>
          <small>{inputPath ? "Use Prepare all folds after editing axis or ROI settings." : "All three symmetry folds will be prepared automatically."}</small>
        </div>
      )}
      {viewerMessage && <div className="viewer-status">{viewerMessage}</div>}
    </div>
  );
}

async function cachedViewerFile(path: string, cache: Map<string, string>) {
  const existing = cache.get(path);
  if (existing) return existing;
  const text = await readViewerFile(path);
  cache.set(path, text);
  return text;
}

async function loadMmcif(plugin: PluginUIContext, text: string, label: string) {
  const data = await plugin.builders.data.rawData({ data: text, label });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, "mmcif");
  const model = await plugin.builders.structure.createModel(trajectory);
  return plugin.builders.structure.createStructure(model);
}

async function addAllRepresentation(
  plugin: PluginUIContext,
  structure: any,
  key: string,
  color: Color,
  alpha: number,
) {
  const component = await plugin.builders.structure.tryCreateComponentStatic(
    structure,
    "all",
    { label: key },
  );
  if (!component) return;
  await plugin.builders.structure.representation.addRepresentation(component, {
    type: "point",
    color: "uniform",
    colorParams: { value: color },
    typeParams: { alpha, pointSize: 2.2 },
  } as any);
}

async function addWeightPreview(plugin: PluginUIContext, structure: any, weight: number) {
  const retained = MS.struct.generator.atomGroups({
    "atom-test": MS.core.rel.gre([
      MS.struct.atomProperty.macromolecular.B_iso_or_equiv(),
      weight,
    ]),
  });
  const removed = MS.struct.generator.atomGroups({
    "atom-test": MS.core.rel.lt([
      MS.struct.atomProperty.macromolecular.B_iso_or_equiv(),
      weight,
    ]),
  });

  const retainedComponent = await plugin.builders.structure.tryCreateComponentFromExpression(
    structure,
    retained,
    "weight-retained",
    { label: "Estimated retained region" },
  );
  if (retainedComponent) {
    await plugin.builders.structure.representation.addRepresentation(retainedComponent, {
      type: "point",
      color: "uniform",
      colorParams: { value: COLORS.retained },
      typeParams: { alpha: 1, pointSize: 2.4 },
    } as any);
  }

  const removedComponent = await plugin.builders.structure.tryCreateComponentFromExpression(
    structure,
    removed,
    "weight-removed",
    { label: "Estimated removed region" },
  );
  if (removedComponent) {
    await plugin.builders.structure.representation.addRepresentation(removedComponent, {
      type: "point",
      color: "uniform",
      colorParams: { value: COLORS.removed },
      typeParams: { alpha: 0.16, pointSize: 1.8 },
    } as any);
  }
}

async function addRetainedRepresentation(
  plugin: PluginUIContext,
  structure: any,
  weight: number,
  key: string,
  color: Color,
) {
  const retained = MS.struct.generator.atomGroups({
    "atom-test": MS.core.rel.gre([
      MS.struct.atomProperty.macromolecular.B_iso_or_equiv(),
      weight,
    ]),
  });
  const component = await plugin.builders.structure.tryCreateComponentFromExpression(
    structure,
    retained,
    key,
    { label: "Estimated sliced region" },
  );
  if (!component) return;
  await plugin.builders.structure.representation.addRepresentation(component, {
    type: "point",
    color: "uniform",
    colorParams: { value: color },
    typeParams: { alpha: 0.95, pointSize: 2.4 },
  } as any);
}

async function addRoiRepresentation(
  plugin: PluginUIContext,
  structure: any,
  asymIds: string[],
  startResid: number,
  endResid: number,
  minimumNormalizedZ?: number,
) {
  const tests: Record<string, any> = {
    "chain-test": MS.core.set.has([
      MS.set(...asymIds),
      MS.struct.atomProperty.macromolecular.label_asym_id(),
    ]),
    "residue-test": MS.core.rel.inRange([
      MS.struct.atomProperty.macromolecular.label_seq_id(),
      startResid,
      endResid,
    ]),
  };
  if (minimumNormalizedZ !== undefined) {
    tests["atom-test"] = MS.core.rel.gre([
      MS.struct.atomProperty.macromolecular.B_iso_or_equiv(),
      minimumNormalizedZ,
    ]);
  }
  const expression = MS.struct.generator.atomGroups(tests);
  const component = await plugin.builders.structure.tryCreateComponentFromExpression(
    structure,
    expression,
    "roi-highlight",
    { label: `ROI ${asymIds.length}-fold copies:${startResid}-${endResid}` },
  );
  if (!component) return;
  await plugin.builders.structure.representation.addRepresentation(component, {
    type: "spacefill",
    color: "uniform",
    colorParams: { value: COLORS.roi },
    typeParams: { alpha: 1, sizeFactor: 0.55 },
  } as any);
}

function getRoiFrames(
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

async function addAsymRepresentation(
  plugin: PluginUIContext,
  structure: any,
  asymId: string,
  key: string,
  color: Color,
  alpha: number,
  sizeFactor: number,
) {
  const expression = MS.struct.generator.atomGroups({
    "chain-test": MS.core.rel.eq([
      MS.struct.atomProperty.macromolecular.label_asym_id(),
      asymId,
    ]),
  });
  const component = await plugin.builders.structure.tryCreateComponentFromExpression(
    structure,
    expression,
    key,
    { label: key },
  );
  if (!component) return;
  await plugin.builders.structure.representation.addRepresentation(component, {
    type: "spacefill",
    color: "uniform",
    colorParams: { value: color },
    typeParams: { alpha, sizeFactor },
  } as any);
}

function getBounds(diagnostics: Record<string, unknown> | undefined): Bounds | null {
  const bounds = diagnostics?.bounds as Bounds | undefined;
  if (!bounds || !Array.isArray(bounds.min) || !Array.isArray(bounds.max)) return null;
  return bounds;
}

function buildMarkerCif(bounds: Bounds, weight: number, axis: boolean, plane: boolean) {
  const [minX, minY, minZ] = bounds.min;
  const [maxX, maxY, maxZ] = bounds.max;
  const centerX = (minX + maxX) / 2;
  const centerY = (minY + maxY) / 2;
  const cutoffZ = minZ + (maxZ - minZ) * weight;
  const rows: string[] = [];
  let id = 1;

  if (axis) {
    for (let index = 0; index <= 64; index += 1) {
      const z = minZ + ((maxZ - minZ) * index) / 64;
      rows.push(markerRow(id++, "AXIS", index + 1, centerX, centerY, z));
    }
  }
  if (plane) {
    const paddingX = (maxX - minX) * 0.08;
    const paddingY = (maxY - minY) * 0.08;
    for (let ix = 0; ix <= 12; ix += 1) {
      for (let iy = 0; iy <= 12; iy += 1) {
        const x = minX + paddingX + ((maxX - minX - 2 * paddingX) * ix) / 12;
        const y = minY + paddingY + ((maxY - minY - 2 * paddingY) * iy) / 12;
        rows.push(markerRow(id++, "PLANE", id, x, y, cutoffZ));
      }
    }
  }

  return `data_CAPSACIN_GUIDES
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
${rows.join("\n")}
#
`;
}

function markerRow(id: number, asymId: string, resid: number, x: number, y: number, z: number) {
  return `HETATM ${id} C C GLY ${asymId} 1 ${resid} C GLY ${asymId} ${resid} ${x.toFixed(3)} ${y.toFixed(3)} ${z.toFixed(3)} 1.00 0.00`;
}
