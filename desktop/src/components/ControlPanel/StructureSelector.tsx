import { useEffect, useState } from "react";
import { useAppStore, BUILTIN_PDBS } from "../../store/useAppStore";
import { preparePreviews } from "../../sidecar/preparePreviews";

export function StructureSelector() {
  const inputPath = useAppStore((s) => s.inputPath);
  const isBuiltin = useAppStore((s) => s.isBuiltin);
  const structureName = useAppStore((s) => s.structureName);
  const setInputPath = useAppStore((s) => s.setInputPath);
  const sidecarStatus = useAppStore((s) => s.sidecarStatus);
  const setError = useAppStore((s) => s.setError);
  const [builtinPdbs, setBuiltinPdbs] = useState(BUILTIN_PDBS);

  useEffect(() => {
    import("../../sidecar/client")
      .then(({ listBuiltinStructures }) => listBuiltinStructures())
      .then((names) => {
        if (names.length > 0) setBuiltinPdbs(names);
      })
      .catch(() => { /* Browser-only preview uses the repository fallback list. */ });
  }, []);

  const handleBuiltinChange = async (name: string) => {
    if (!name) return;
    try {
      const { getStructurePath } = await import("../../sidecar/client");
      let path: string;
      try {
        path = await getStructurePath(name);
      } catch {
        // Browser-only development mode does not expose Tauri's resource resolver.
        path = name;
      }
      setInputPath(path, true, name);
      await preparePreviews();
    } catch (error: any) {
      setError(error?.message ?? String(error));
    }
  };

  const handleOpenFile = async () => {
    try {
      const { openPdbDialog } = await import("../../sidecar/client");
      const path = await openPdbDialog();
      if (path) {
        const name = path.split("/").pop()?.replace(".pdb", "") ?? "imported";
        setInputPath(path, false, name);
        await preparePreviews();
      }
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    }
  };

  const busy = sidecarStatus === "running";

  return (
    <div className="panel-section">
      <h3>Structure</h3>
      <select
        value={isBuiltin ? structureName : ""}
        onChange={(e) => handleBuiltinChange(e.target.value)}
        disabled={busy}
        style={{ marginBottom: "8px" }}
      >
        <option value="">— Select built-in —</option>
        {builtinPdbs.map((name) => (
          <option key={name} value={name}>{name}.pdb</option>
        ))}
      </select>
      <button
        className="btn btn-outline"
        onClick={handleOpenFile}
        disabled={busy}
        style={{ width: "100%", fontSize: "12px" }}
      >
        Open Local PDB…
      </button>
      {inputPath && (
        <div style={{ marginTop: "8px", fontSize: "11px", color: "var(--text-secondary)" }}>
          {isBuiltin ? `📦 ${structureName}.pdb` : `📂 ${inputPath.split("/").pop()}`}
        </div>
      )}
    </div>
  );
}
