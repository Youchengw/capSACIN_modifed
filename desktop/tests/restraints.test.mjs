import { beforeEach, test } from "node:test";
import assert from "node:assert/strict";
import { useAppStore } from "../src/store/useAppStore.ts";
import { generateRestraints } from "../src/sidecar/generateRestraints.ts";

const state = () => useAppStore.getState();
const slice = () => ({ output_pdb_path: "/tmp/sliced.pdb", sliced_cif_path: "/tmp/sliced.cif" });
const result = () => ({ input_path: "/tmp/sliced.pdb", archive_path: "/tmp/restraints.zip" });
beforeEach(() => {
  state().reset();
  state().setInputPath("/tmp/input.pdb", false, "input");
  state().setSliceResult(slice());
});

test("generates from the final slice with independent strength and selects the region view", async () => {
  state().setWeight(0.7);
  state().setSliceResult(slice());
  let calls = 0;
  await generateRestraints({
    onProgress: async () => () => {},
    generateRestraints: async (params) => {
      calls++;
      assert.deepEqual(params, { input_path: "/tmp/sliced.pdb", full_until: 15, free_from: 30, strength: 1 });
      return result();
    },
  });
  assert.equal(state().displayMode, "restraints");
  assert.ok(state().restraintResult);
  state().setDisplayMode("sliced");
  state().setDisplayMode("restraints");
  assert.equal(calls, 1);
});

test("geometry or restraint settings invalidate old downloads; display-only edits preserve them", () => {
  for (const edit of [
    () => state().setWeight(0.8), () => state().setSymmetry(3),
    () => state().setAxisIndex(1), () => state().setRoiEnabled(true),
    () => state().setInputPath("/tmp/new.pdb", false, "new"),
    () => state().setSliceResult(slice()), () => state().setRestraintSettings({ strength: 2 }),
  ]) {
    state().setRestraintResult(result());
    edit();
    assert.equal(state().restraintResult, null);
  }
  state().setRoiEnabled(false);
  const saved = result();
  state().setRestraintResult(saved);
  state().setRoiStartResid(42);
  state().toggleAxis();
  assert.equal(state().restraintResult, saved);
});

test("a response for an invalidated slice cannot restore a stale archive", async () => {
  for (const edit of [() => state().setWeight(0.9),
    () => state().setInputPath("/tmp/new.pdb", false, "new"),
    () => state().setRestraintSettings({ free_from: 40 })]) {
    state().setSliceResult(slice());
    let finish;
    const pending = generateRestraints({ onProgress: async () => () => {},
      generateRestraints: () => new Promise((resolve) => { finish = resolve; }) });
    await new Promise((resolve) => setImmediate(resolve));
    edit();
    finish(result());
    await pending;
    assert.equal(state().restraintResult, null);
    assert.equal(state().sidecarStatus, "idle");
  }
});

test("failed or cancelled generation leaves no downloadable result and allows retry", async () => {
  for (const message of ["Operation cancelled", "No usable atoms"]) {
    state().setRestraintResult(result());
    await generateRestraints({ onProgress: async () => () => {},
      generateRestraints: async () => { throw new Error(message); } });
    assert.equal(state().restraintResult, null);
    assert.equal(state().sidecarStatus, message.includes("cancelled") ? "cancelled" : "error");
  }
});
