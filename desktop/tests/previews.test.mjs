import { beforeEach, test } from "node:test";
import assert from "node:assert/strict";
import { useAppStore } from "../src/store/useAppStore.ts";
import { preparePreviews } from "../src/sidecar/preparePreviews.ts";

const store = () => useAppStore.getState();
const preview = (fold) => ({
  axis_candidates: [], selected_axis_index: 0, selected_axis: [0, 0, 1],
  aligned_cif_path: `/tmp/prepared-${fold}.cif`, workspace_path: `/tmp/fold-${fold}`,
  diagnostics: {}, warnings: [],
});
const batch = () => ({
  inspection: { input_path: "/tmp/capsid.pdb", n_models: 60, chains: [],
    n_chains_per_frame: 1, total_atoms: 100, total_residues: 10, valid_residue_ranges: {} },
  previews: { 2: preview(2), 3: preview(3), 5: preview(5) }, errors: {},
});
function client(result = batch()) {
  const calls = [];
  return {
    calls,
    onProgress: async () => () => {},
    prepareAllPreviews: async (params) => { calls.push(params); return result; },
    preparePreview: async (params) => { calls.push(params); return preview(params.symmetry); },
  };
}
beforeEach(() => {
  store().reset();
  store().setInputPath("/tmp/capsid.pdb", false, "capsid");
});

test("opening prepares all folds once; 5 -> 3 -> 2 -> 5 and weight edits need no prepare", async () => {
  const api = client();
  await preparePreviews(api);
  const five = store().previewResult;
  for (const fold of [3, 2, 5]) {
    store().setSymmetry(fold);
    assert.equal(store().previewResult.aligned_cif_path, `/tmp/prepared-${fold}.cif`);
  }
  store().setWeight(0.8);
  assert.equal(store().previewResult, five);
  assert.equal(api.calls.length, 1);
  assert.equal(store().sidecarStatus, "idle");
});

test("candidate changes invalidate only that fold and remember its selection", async () => {
  const api = client();
  await preparePreviews(api);
  store().setAxisIndex(2);
  assert.equal(store().previewResult, null);
  store().setSymmetry(3);
  assert.equal(store().axisIndex, 0);
  assert.ok(store().previewResult);
  store().setSymmetry(5);
  assert.equal(store().axisIndex, 2);
  assert.equal(store().previewResult, null);
  await preparePreviews(api);
  assert.deepEqual(api.calls[1].axis_indices, { 2: 0, 3: 0, 5: 2 });
});

test("ROI changes and another PDB cannot reuse old fold previews", async () => {
  await preparePreviews(client());
  store().setRoiEnabled(true);
  store().setRoiSelection("protein and resid 42:50");
  store().setSymmetry(3);
  assert.equal(store().previewResult, null);
  assert.deepEqual(store().foldPreviews, {});
  await preparePreviews(client());
  store().setInputPath("/tmp/other.pdb", false, "other");
  assert.deepEqual(store().foldPreviews, {});
  assert.equal(store().previewResult, null);
});

test("editing a display-only ROI preserves all prepared folds", async () => {
  await preparePreviews(client());
  const previews = store().foldPreviews;
  store().setRoiStartResid(42);
  store().setRoiEndResid(50);
  store().setRoiFrame(1);
  assert.equal(store().foldPreviews, previews);
  store().setSymmetry(3);
  assert.equal(store().previewResult, previews[3]);
});

test("a failed fold stays unavailable while successful folds can be selected", async () => {
  const result = batch();
  delete result.previews[3];
  result.errors[3] = "No usable axis";
  await preparePreviews(client(result));
  assert.ok(store().previewResult);
  store().setSymmetry(3);
  assert.equal(store().previewResult, null);
  assert.equal(store().error, "No usable axis");
  store().setSymmetry(2);
  assert.ok(store().previewResult);
  assert.equal(store().error, null);
});

test("a late response for the previous PDB cannot replace the new structure", async () => {
  let finish;
  const api = client();
  api.prepareAllPreviews = () => new Promise((resolve) => { finish = resolve; });
  const pending = preparePreviews(api);
  await new Promise((resolve) => setImmediate(resolve));
  store().setInputPath("/tmp/other.pdb", false, "other");
  finish(batch());
  await pending;
  assert.equal(store().inputPath, "/tmp/other.pdb");
  assert.equal(store().previewResult, null);
  assert.equal(store().inspectResult, null);
});

test("cancellation permits retry and manual mode prepares only the selected fold", async () => {
  const api = client();
  api.prepareAllPreviews = async () => { throw new Error("Operation cancelled"); };
  await preparePreviews(api);
  assert.equal(store().sidecarStatus, "cancelled");
  assert.equal(store().previewResult, null);
  store().setAxisMode("manual");
  store().setRefIndices("1,2,3");
  const manualApi = client();
  await preparePreviews(manualApi);
  assert.equal(manualApi.calls.length, 1);
  assert.equal(manualApi.calls[0].auto, false);
  assert.deepEqual(Object.keys(store().foldPreviews), ["5"]);
});
