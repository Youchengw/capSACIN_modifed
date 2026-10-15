import { beforeEach, test } from "node:test";
import assert from "node:assert/strict";
import { useAppStore } from "../src/store/useAppStore.ts";
import { applyChainSelection } from "../src/sidecar/applyChainSelection.ts";
import { chainSelectionDirty, weightChainIndices } from "../src/utils/chains.ts";

const state = () => useAppStore.getState();
const candidate = (source_index, z_mean, selected) => ({
  id: `M${source_index}_X`, source_index, frame: source_index, original_chain_id: "X",
  z_min: z_mean - 5, z_mean, z_max: z_mean + 5, atom_count: 2, residue_count: 1,
  selected, selected_by_weight: selected,
});
const slice = (selected = [1]) => ({
  output_pdb_path: "/tmp/sliced.pdb", sliced_cif_path: "/tmp/sliced.cif",
  workspace_path: "/tmp/slice", chain_source_path: "/tmp/chains.pkl",
  chain_candidates: [candidate(0, 10, selected.includes(0)), candidate(1, 20, selected.includes(1))],
  statistics: { input_atoms: 4, output_atoms: selected.length * 2, input_residues: 2,
    output_residues: selected.length, input_chains: 2, output_chains: selected.length,
    retained_percentage: selected.length * 50 }, axis: [0, 0, 1], diagnostics: {}, warnings: [],
});

beforeEach(() => {
  state().reset();
  state().setInputPath("/tmp/input.pdb", false, "input");
  state().setSliceResult(slice());
});

test("chain edits track the weight result and invalidate restraints", () => {
  assert.deepEqual(state().selectedChainIndices, [1]);
  assert.deepEqual(weightChainIndices(state().sliceResult), [1]);
  assert.equal(chainSelectionDirty(state().sliceResult, state().selectedChainIndices), false);
  state().setRestraintResult({ archive_path: "/tmp/restraints.zip" });
  state().setChainSelected(0, true);
  assert.deepEqual(state().selectedChainIndices.sort(), [0, 1]);
  assert.equal(state().restraintResult, null);
  assert.equal(chainSelectionDirty(state().sliceResult, state().selectedChainIndices), true);
});

test("applying selection replaces the exported model and preserves weight markers", async () => {
  state().setChainSelected(0, true);
  state().setChainSelected(1, false);
  let params;
  await applyChainSelection({
    onProgress: async () => () => {},
    applyChainSelection: async (value) => { params = value; return slice([0]); },
  });
  assert.deepEqual(params, { chain_source_path: "/tmp/chains.pkl", selected_chain_indices: [0] });
  assert.deepEqual(state().selectedChainIndices, [0]);
  assert.equal(state().sliceResult.statistics.output_chains, 1);
  assert.equal(chainSelectionDirty(state().sliceResult, state().selectedChainIndices), false);
});
