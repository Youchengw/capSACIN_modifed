#!/usr/bin/env python3
"""
Inspect residue ID ranges in CapSACIN-style PDB files.

CapSACIN inputs often store the 60 capsid copies as MODEL records, where each
MODEL contains one subunit/asymmetric-unit-like copy.  This helper reports the
residue ID ranges per MODEL and chain so users can choose valid
``--roi-selection`` values for ``sliceCapsid.py``.
"""

import argparse
from collections import defaultdict
from pathlib import Path


def resolve_pdb_path(pdb_arg):
    path = Path(pdb_arg)
    if path.exists():
        return path

    if path.suffix.lower() != ".pdb":
        candidates = [
            Path("input") / f"{pdb_arg}.pdb",
            Path(__file__).resolve().parent / "input" / f"{pdb_arg}.pdb",
        ]
    else:
        candidates = [
            Path("input") / path.name,
            Path(__file__).resolve().parent / "input" / path.name,
        ]

    for candidate in candidates:
        if candidate.exists():
            return candidate

    tried = ", ".join(str(c) for c in candidates)
    raise FileNotFoundError(f"Could not find PDB {pdb_arg!r}. Tried: {tried}")


def parse_residue_ranges(pdb_path, include_hetatm=False):
    records = ("ATOM", "HETATM") if include_hetatm else ("ATOM",)
    model_number = None
    implicit_model_started = False
    ranges = defaultdict(lambda: defaultdict(set))

    with pdb_path.open() as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                try:
                    model_number = int(line[10:14])
                except ValueError:
                    model_number = len(ranges) + 1
                continue
            if record == "ENDMDL":
                model_number = None
                continue
            if record not in records:
                continue

            if model_number is None:
                model_number = 1
                implicit_model_started = True

            chain_id = line[21].strip() or "-"
            resid_text = line[22:26].strip()
            insertion_code = line[26].strip()
            if not resid_text:
                continue
            try:
                resid = int(resid_text)
            except ValueError:
                continue

            ranges[model_number][chain_id].add((resid, insertion_code))

            if implicit_model_started:
                model_number = 1

    return ranges


def summarize_range(residue_keys):
    numeric_resids = [resid for resid, _ in residue_keys]
    insertion_codes = sorted({icode for _, icode in residue_keys if icode})
    return {
        "min": min(numeric_resids),
        "max": max(numeric_resids),
        "count": len(set(residue_keys)),
        "insertions": insertion_codes,
    }


def format_insertions(insertions):
    if not insertions:
        return "-"
    return ",".join(insertions)


def print_summary(pdb_path, ranges, per_model_limit=None):
    print(f"PDB: {pdb_path}")
    print(f"Models: {len(ranges)}")
    print("")

    global_by_chain = defaultdict(set)
    for chain_map in ranges.values():
        for chain_id, residue_keys in chain_map.items():
            global_by_chain[chain_id].update(residue_keys)

    print("Global residue ranges by chain:")
    print("chain  min_resid  max_resid  unique_resids  insertion_codes")
    for chain_id in sorted(global_by_chain):
        summary = summarize_range(global_by_chain[chain_id])
        print(
            f"{chain_id:>5}  {summary['min']:>9}  {summary['max']:>9}  "
            f"{summary['count']:>13}  {format_insertions(summary['insertions'])}"
        )

    print("")
    print("Per-MODEL residue ranges:")
    print("model  roi_frame  chain  min_resid  max_resid  unique_resids")

    model_numbers = sorted(ranges)
    if per_model_limit is not None:
        model_numbers = model_numbers[:per_model_limit]

    for model_number in model_numbers:
        roi_frame = model_number - 1
        for chain_id in sorted(ranges[model_number]):
            summary = summarize_range(ranges[model_number][chain_id])
            print(
                f"{model_number:>5}  {roi_frame:>9}  {chain_id:>5}  "
                f"{summary['min']:>9}  {summary['max']:>9}  {summary['count']:>13}"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Report residue ID ranges in a PDB file."
    )
    parser.add_argument(
        "pdb",
        help="PDB path or input prefix, for example '9jjh' or 'input/9jjh.pdb'.",
    )
    parser.add_argument(
        "--include-hetatm",
        action="store_true",
        help="Include HETATM records in addition to ATOM records.",
    )
    parser.add_argument(
        "--per-model-limit",
        type=int,
        default=5,
        help="Number of MODEL entries to print in detail. Use 0 for all.",
    )
    args = parser.parse_args()

    pdb_path = resolve_pdb_path(args.pdb)
    ranges = parse_residue_ranges(pdb_path, include_hetatm=args.include_hetatm)
    if not ranges:
        raise ValueError(f"No residue records found in {pdb_path}")

    per_model_limit = args.per_model_limit
    if per_model_limit == 0:
        per_model_limit = None
    print_summary(pdb_path, ranges, per_model_limit=per_model_limit)


if __name__ == "__main__":
    main()
