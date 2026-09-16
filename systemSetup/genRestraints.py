#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
from capsacin.restraints import generate_restraint_files


def main(args):
    title, weight, symmetry = args.pdb, args.weight, args.symmetry
    suffix = f"-sym{symmetry}" if symmetry is not None else ""
    result = generate_restraint_files(
        f"./output/{title}-sliced{suffix}-w{weight}.pdb",
        f"output/{title}-w{weight}-posre", strength=weight, prefix=title,
    )
    restrained = result["regions"] != 3
    plt.scatter(result["z"][restrained], result["force_constants"][restrained])
    plt.xlabel("z (Å)")
    plt.ylabel("Force constant (kJ mol⁻¹ nm⁻²)")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=str, required=True, help="PDB filename prefix used for sliceCapsid.py")
    parser.add_argument("--weight", type=float, default=1.0, help="Restraint scaling weight (default 1.0)")
    parser.add_argument("--symmetry", "--sym", type=int, choices=(2, 3, 5), default=None,
                        help="Symmetry used for the sliced PDB name. If omitted, reads the legacy name.")
    args = parser.parse_args()
    main(args)
