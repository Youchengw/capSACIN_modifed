import { Mesh } from "molstar/lib/mol-geo/geometry/mesh/mesh";
import { MeshBuilder } from "molstar/lib/mol-geo/geometry/mesh/mesh-builder";
import { Vec3 } from "molstar/lib/mol-math/linear-algebra";
import { Shape } from "molstar/lib/mol-model/shape";
import type { Structure } from "molstar/lib/mol-model/structure";
import { PluginStateObject as SO, PluginStateTransform } from "molstar/lib/mol-plugin-state/objects";
import { ShapeRepresentation3D } from "molstar/lib/mol-plugin-state/transforms/representation";
import type { PluginUIContext } from "molstar/lib/mol-plugin-ui/context";
import { Color } from "molstar/lib/mol-util/color";
import { ParamDefinition as PD } from "molstar/lib/mol-util/param-definition";
import type { RestraintSettings } from "../../sidecar/types";

const RestraintPlaneShape = PluginStateTransform.BuiltIn({
  name: "capsacin-restraint-boundary-planes",
  display: { name: "Restraint boundaries" },
  from: SO.Root,
  to: SO.Shape.Provider,
  params: {
    min: PD.Vec3(Vec3()),
    max: PD.Vec3(Vec3()),
    fullUntil: PD.Numeric(15),
    freeFrom: PD.Numeric(30),
  },
})({
  apply({ params }) {
    return new SO.Shape.Provider({
      label: "Restraint boundaries",
      data: params,
      params: Mesh.Params,
      geometryUtils: Mesh.Utils,
      getShape: (_, data) => {
        const mesh = MeshBuilder.createState(12, 12);
        const padding = Math.max(3, Math.max(data.max[0] - data.min[0], data.max[1] - data.min[1]) * 0.04);
        const x0 = data.min[0] - padding;
        const y0 = data.min[1] - padding;
        const x1 = data.max[0] + padding;
        const y1 = data.max[1] + padding;
        const heights = [data.fullUntil, data.freeFrom];
        heights.forEach((z, group) => {
          mesh.currentGroup = group;
          const a = Vec3.create(x0, y0, z);
          const b = Vec3.create(x1, y0, z);
          const c = Vec3.create(x1, y1, z);
          const d = Vec3.create(x0, y1, z);
          MeshBuilder.addTriangle(mesh, a, b, c);
          MeshBuilder.addTriangle(mesh, a, c, d);
        });
        return Shape.create(
          "Restraint boundaries", data, MeshBuilder.getMesh(mesh),
          (group) => Color(group === 0 ? 0xfb7185 : 0xfbbf24),
          () => 1,
          (group) => `${group === 0 ? "Full / transition" : "Transition / free"} · z = ${heights[group]} Å`,
        );
      },
    }, { label: "Restraint boundaries" });
  },
});

export async function addRestraintPlanes(
  plugin: PluginUIContext,
  structure: Structure,
  settings: RestraintSettings,
) {
  // Use the final sliced coordinates, whose z origin matches the exported ITPs.
  const bounds = structure.boundary.box;
  const representation = await plugin.build().toRoot().apply(RestraintPlaneShape, {
    min: Vec3.clone(bounds.min), max: Vec3.clone(bounds.max),
    fullUntil: settings.full_until, freeFrom: settings.free_from,
  }).apply(ShapeRepresentation3D, {
    alpha: 0.2, doubleSided: true, transparentBackfaces: "on", ignoreLight: true,
  }).commit();
  // The guides should not intercept atom picking beneath their surfaces.
  representation.obj?.data.repr.setState({ pickable: false });
  return representation.ref;
}
