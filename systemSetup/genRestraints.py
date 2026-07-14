#!/usr/bin/env python3
import argparse
import numpy as np
from numpy import savetxt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import MDAnalysis as md
from tqdm import tqdm
import pandas as pd
import copy
import os

def main(args):
    title = args.pdb
    weight = args.weight  # For consistency, but note the original code uses z ranges 15-30 as hardcoded
    symmetry = args.symmetry

    if symmetry is None:
        gro = f"./output/{title}-sliced-w{weight}.pdb"
    else:
        gro = f"./output/{title}-sliced-sym{symmetry}-w{weight}.pdb"
    u = md.Universe(gro)

    chainID = list(np.unique(u.select_atoms("protein").chainIDs))
    numChains = len(chainID)

    chains, names, indices, positions = [], [], [], []
    for i in range(numChains):
        currChain = u.select_atoms(f"chainid {chainID[i]}")
        chains.append(np.repeat(chainID[i], len(currChain)))
        names.append(currChain.names)
        indices.append(np.arange(1,len(currChain)+1,1))
        positions.append(currChain.positions)

    chainStack = np.hstack(chains)
    nameStack = np.hstack(names)
    idxStack = np.hstack(indices)
    positions = np.vstack(positions)
    z = np.round(positions[:,2], 3)

    allChains = np.vstack((idxStack, chainStack, nameStack, z)).T
    df = pd.DataFrame(allChains, columns=["idx","chain","name","z"])

    # Normalize and apply sigmoidal function to get restraint weights
    def norm(x):
        num = (x - x.min()) * 2.0
        denom = x.max() - x.min()
        return num/denom - 1.0

    def sig(x, k=0.1):
        s = -1/(1+np.exp(-x/k)) + 1  # inverted sigmoid
        return s

    # Sigmoidal weight for everything between 15 and 30 angstroms
    selection = df[(df['z'] < 30) & (df['z'] >= 15)]
    z_norm = np.asarray(norm(selection.z), dtype='float')
    s = sig(z_norm)
    plt.scatter(z_norm, s)
    plt.show()
    s = pd.DataFrame(s, columns=['s'])
    selection = selection.reset_index(drop=True)
    selectionSigmoid = pd.concat((selection, s), axis=1)

    # Everything between 0 < z < 15 gets full weight
    selection2 = df[df['z'] < 15]
    s2 = np.linspace(1,1,len(selection2))
    s2 = pd.DataFrame(s2, columns=['s'])
    selection2 = selection2.reset_index(drop=True)
    selection2Sigmoid = pd.concat((selection2, s2), axis=1)

    # Concatenate
    fullWeights = pd.concat((selectionSigmoid, selection2Sigmoid), axis=0)
    K = fullWeights.s * 1000 * weight  # Multiply by weight from CLI
    K = K.rename("K")
    fullWeights = pd.concat((fullWeights, K), axis=1)
    plt.clf()
    plt.scatter(fullWeights.z, fullWeights.K)

    # Iterate through chains to create individual POSRES files
    header="[ position_restraints ] \n; atom  type      fx      fy      fz"
    for chain in chainID:
        group = fullWeights[fullWeights['chain'] == chain]
        idx = group.idx
        K_chain = np.round(group.K, 3)
        typeList = np.linspace(1,1,len(idx))
        output = np.vstack((idx, typeList, K_chain, K_chain, K_chain)).T
        output = output[output[:,0].argsort()]

        output_dir = f'output/{title}-w{weight}-posre/'
        os.makedirs(output_dir, exist_ok=True)

        savetxt(f'{output_dir}{title}-posre-scaled_Protein_chain_{chain}.itp', 
                output, fmt='%s', header=header, comments='', delimiter=' ')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=str, required=True, help="PDB filename prefix used for sliceCapsid.py")
    parser.add_argument("--weight", type=float, default=1.0, help="Restraint scaling weight (default 1.0)")
    parser.add_argument("--symmetry", "--sym", type=int, choices=(2, 3, 5), default=None,
                        help="Symmetry used for the sliced PDB name. If omitted, reads the legacy name.")
    args = parser.parse_args()
    main(args)
