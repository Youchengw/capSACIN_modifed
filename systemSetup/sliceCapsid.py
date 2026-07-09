#!/usr/bin/env python3
import argparse
import warnings
import numpy as np
from numpy import savetxt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import MDAnalysis as md
from tqdm import tqdm
import pandas as pd
from capsacin import definePlane, formatPDB, createDictionary, findSymmetryAxes
import copy
import os

if not os.path.exists("output"):
    os.makedirs("output")

# Suppress pandas SettingWithCopyWarning safely
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)

def split_dataframe(df, chunk_size): 
    chunks = list()
    num_chunks = len(df) // chunk_size + 1
    for i in range(num_chunks):
        chunks.append(df[i*chunk_size:(i+1)*chunk_size])
    return chunks

def rotationMatrix(vec1, vec2):
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def print_axis_table(rows, has_roi=False):
    if len(rows) == 0:
        print("[axes] No axes found.")
        return

    if has_roi:
        print("[axes] rank axis-index roi-score roi-angle-deg roi-line-dist symmetry-score axis")
        for rank, row in enumerate(rows):
            axis = row["axis"]
            print(
                "[axes] "
                f"{rank:>4} {row['axis_index']:>10} "
                f"{row['roi_score']:.5f} {row['roi_angle_deg']:.3f} "
                f"{row['roi_line_distance']:.3f} {row['score']:.5f} "
                f"[{axis[0]:.5f}, {axis[1]:.5f}, {axis[2]:.5f}]"
            )
    else:
        print("[axes] rank axis-index symmetry-score ref-index ref-frame chain axis")
        for rank, row in enumerate(rows):
            axis = row["axis"]
            print(
                "[axes] "
                f"{rank:>4} {row['axis_index']:>10} {row['score']:.5f} "
                f"{row['ref_index']:>9} {row['ref_frame']:>9} "
                f"{row['chain_id']} "
                f"[{axis[0]:.5f}, {axis[1]:.5f}, {axis[2]:.5f}]"
            )

def main(args):
    # CLI options
    output_title = args.pdb
    symmetry = args.symmetry
    Plot = args.plot
    indicesVMD = args.refindex
    weight = args.weight
    brokenResCheck = True
    rotationCheck = True
    auto_mode = args.auto
    axis_index = args.axis_index
    roi_selection = args.roi_selection
    roi_frame = args.roi_frame

    gro=f"input/{output_title}.pdb"
    traj=f"input/{output_title}.pdb"
    u=md.Universe(gro, traj)

    if args.list_axes:
        axis_rows = findSymmetryAxes.list_axes(
            u,
            symmetry,
            roi_selection=roi_selection,
            roi_frame=roi_frame,
        )
        print_axis_table(axis_rows, has_roi=bool(roi_selection))
        if not auto_mode:
            return

    if roi_selection and not auto_mode:
        print("[roi] --roi-selection is only used with --auto or --list-axes.")

    # --- Auto-detect symmetry axis and reference atoms ---
    auto_ref_indices = None  # [pointA_idx, pointB_idx, pointC_idx]
    if auto_mode:
        auto_ref_indices, axis_dir, auto_frames = findSymmetryAxes.auto_detect_reference(
            u,
            symmetry,
            axis_index=axis_index,
            roi_selection=roi_selection,
            roi_frame=roi_frame,
        )
        indicesVMD = [auto_ref_indices[0]]  # pointA for backward compat
        if roi_selection:
            print(f"[auto] ROI selection: {roi_selection!r} on frame {roi_frame}")
            print(f"[auto] Using ROI-ranked {symmetry}-fold axis at rank {axis_index}")
        print(f"[auto] Detected {symmetry}-fold axis: {axis_dir}")
        print(f"[auto] Reference atom indices: {auto_ref_indices}")
        print(f"[auto] Reference frames: {auto_frames}")
    # ----------------------------------------------------

    origChains = np.unique(u.select_atoms("protein").chainIDs)
    nMonomers = len(origChains)

    # Reference position
    refPos = u.select_atoms(f"protein and index {indicesVMD[0]}")
    pointA = refPos.positions[0]

    # In auto mode, get pointB and pointC directly from the detected monomers.
    # This bypasses the nearest-neighbour distance heuristic, which can fail
    # for certain atom choices.
    if auto_mode and auto_ref_indices is not None:
        # pointB from the second detected monomer
        u.trajectory[auto_frames[1]]
        refPosB_md = u.select_atoms(f"protein and index {auto_ref_indices[1]}")
        pointB = refPosB_md.positions[0].copy()
        # pointC from the third detected monomer
        u.trajectory[auto_frames[2]]
        refPosC_md = u.select_atoms(f"protein and index {auto_ref_indices[2]}")
        pointC = refPosC_md.positions[0].copy()
        # Return to first monomer for pointA
        u.trajectory[auto_frames[0]]
        refPosA_md = u.select_atoms(f"protein and index {auto_ref_indices[0]}")
        pointA = refPosA_md.positions[0].copy()
        # Apply the same z-perturbation as the manual path for consistency
        pointB[2] += 0.1
        pointC[2] -= 0.1

    # Dictionary for broken residue check
    aminoAcidChecks = createDictionary.createDictionary(u)

    # Build virus coords
    allCoords=[]
    for ts in tqdm(u.trajectory):
        for chain in origChains:
            virusCoords = u.select_atoms(f"protein and chainid {chain}").positions
            allCoords.append(virusCoords)
    builtVirus = np.vstack(allCoords)
    com = (np.mean(builtVirus[:,0]), np.mean(builtVirus[:,1]), np.mean(builtVirus[:,2]))
    centerDist = np.sqrt((builtVirus[:,0]-com[0])**2 + (builtVirus[:,1]-com[1])**2 + (builtVirus[:,2]-com[2])**2)

    if Plot:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(builtVirus[0:-1:50,0], builtVirus[0:-1:50,1], builtVirus[0:-1:50,2], c=centerDist[0:-1:50])
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        plt.tight_layout()
        plt.show()

    # Get PDB information
    j=0
    resids, resnames, bfactors, names, types, occupancies, indices, atoms = [], [], [], [], [], [], [], []
    lenChain=[]
    for chain in origChains:
        virus = u.select_atoms(f"protein and chainid {chain}")
        resids.append(virus.resids)
        resnames.append(virus.resnames)
        bfactors.append(virus.tempfactors)
        names.append(virus.names)
        types.append(virus.types)
        occupancies.append(virus.occupancies)
        indices.append(np.arange(1,len(allCoords[j])+1, 1))
        lenChain.append(len(np.unique(resids[j])))
        atom = []
        for i in range(len(allCoords[j])):
            atom.append("ATOM")
        atoms.append(atom)
        j+=1

    N = nMonomers * (u.trajectory.n_frames)
    chainID = [chr(65 + i) for i in range(N)]

    # Stack coords
    chains=[]
    k=0
    for i in range(len(allCoords)):
        x = allCoords[i][:,0]
        y = allCoords[i][:,1]
        z = allCoords[i][:,2]
        
        chainList = []
        for j in range(len(allCoords[k])):
            chainList.append(chainID[i])
            
        index = indices[k] + (i*(len(allCoords[k])))
        currChain = np.vstack((atoms[k], index, names[k], resnames[k], chainList, 
                               resids[k], x,y,z, occupancies[k], bfactors[k], types[k])).T
        chains.append(currChain)
        k+=1
        if k == len(origChains):
            k=0

    chains = np.vstack(chains)
    chains = pd.DataFrame(chains, columns=["atom","idx","name","resname","chain",
                                           "resids", "x","y","z", "occ","b","type"])

    xChain = chains.x
    yChain = chains.y
    zChain = chains.z

    coords = np.vstack((xChain, yChain, zChain)).T

    # Symmetric reference points (manual mode only)
    if not (auto_mode and auto_ref_indices is not None):
        candidatePoints = chains[(chains.resids == refPos.resids[0]) & (chains.resname == refPos.resnames[0]) & (chains.name == refPos.names[0])]
        candidateCoords = np.vstack((candidatePoints.x, candidatePoints.y, candidatePoints.z)).T
        from MDAnalysis.analysis import contacts
        distMatrix = contacts.distance_array(pointA, candidateCoords)
        distMatrix = np.vstack((candidatePoints.idx, distMatrix)).T
        distMatrixS = np.argsort(distMatrix[:,1])

        if symmetry == 3:
            k1, k2 = 1, 2
        elif symmetry == 5:
            k1, k2 = 1, 3
        elif symmetry == 2:
            k1, k2 = 0, 1

        refPosB = candidatePoints[(candidatePoints.idx == candidatePoints.idx.values[distMatrixS[k1]])]
        refPosC = candidatePoints[(candidatePoints.idx == candidatePoints.idx.values[distMatrixS[k2]])]
        pointB = np.reshape(np.array(np.vstack((refPosB.x, refPosB.y, refPosB.z+0.1)).T, dtype=float).T, (3,))
        pointC = np.reshape(np.array(np.vstack((refPosC.x, refPosC.y, refPosC.z-0.1)).T, dtype=float).T, (3,))

    # Alignment
    points = pd.DataFrame(np.vstack((pointA, pointB, pointC)), columns=['x', 'y', 'z'])
    x,y,z = points.x, points.y, points.z
    selectPoints = np.vstack((x,y,z)).T
    com = (np.mean(x), np.mean(y), np.mean(z))
    coords = coords - com
    selectPoints = selectPoints - com
    x, y, z = selectPoints[:,0], selectPoints[:,1], selectPoints[:,2]

    # In auto mode, use the detected axis directly as the alignment normal.
    # This is more robust than the 3-point plane method since the axis was
    # found by a global symmetry search.
    if auto_mode and auto_ref_indices is not None:
        normalVector = axis_dir
        print(f"[auto] Using detected axis as alignment normal (bypassing 3-point method)")
    else:
        normalVector = definePlane.definePlane(x,y,z)
    zAxis = [0,0,1]

    if rotationCheck:
        A = rotationMatrix(normalVector, zAxis).T
        newCoords = np.asarray((coords @ A),dtype=float)
        newSelectPoints = np.asarray((selectPoints @ A),dtype=float)
    else:
        newCoords = coords
        newSelectPoints = selectPoints

    # Plot after alignment
    if Plot:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        normalVectorLine = np.asarray((normalVector @ A),dtype=float)
        normalVectorLine = np.vstack(([0,0,0], normalVectorLine))
        oldVectorLine = np.vstack(([0,0,0], normalVector))
        zLine = np.vstack(([0,0,0], zAxis))
        ax.plot(oldVectorLine[:,0], oldVectorLine[:,1], oldVectorLine[:,2], color='firebrick')
        ax.plot(normalVectorLine[:,0], normalVectorLine[:,1], normalVectorLine[:,2], color='indigo')
        ax.plot(zLine[:,0], zLine[:,1], zLine[:,2], color='darkgreen', lw=4)
        ax.scatter(newSelectPoints[:,0], newSelectPoints[:,1], newSelectPoints[:,2], color="indigo", s=40)
        ax.scatter(selectPoints[:,0], selectPoints[:,1], selectPoints[:,2], color="firebrick", s=40)
        lim=10
        ax.set_xlim(-1*lim,lim)
        ax.set_ylim(-1*lim,lim)
        ax.set_zlim(-1*lim,lim)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        plt.tight_layout()
        plt.show()

    # Center coordinates and prepare for slicing
    newCoords = newCoords - np.mean(newCoords, axis=0)
    newCoords = newCoords - [0,0,np.min(newCoords[:,2])]
    newCoords[:,0] -= np.min(newCoords[:,0])
    newCoords[:,1] -= np.min(newCoords[:,1])
    newCoords = pd.DataFrame(newCoords, columns=["x","y","z"])

    # Update chains dataframe
    chains['x'] = newCoords['x'].apply('{0:.3f}'.format)
    chains['y'] = newCoords['y'].apply('{0:.3f}'.format)
    chains['z'] = newCoords['z'].apply('{0:.3f}'.format)
    chains['b'] = chains['b'].astype('float').apply('{0:.2f}'.format)
    chains['occ'] = chains['occ'].astype('float').apply('{0:.2f}'.format)
    chains = chains.astype('string')
    chains = formatPDB.formatPDB(chains)

    # Slice
    cutoff = np.max(chains.z.astype(float))*weight
    originalLen = len(chains)
    slicedChains = np.unique(chains.chain[chains["z"].astype(float) >= cutoff])
    originalLengths=[]
    for x in slicedChains:
        originalLengths.append(len(np.unique(chains.resids[chains.chain == x])))
    slicedCoords = chains[chains["z"].astype(float) >= cutoff]
    slicedCoords['z'] = slicedCoords['z'].astype(float) - np.min(slicedCoords['z'].astype(float))

    # Broken chains / residues cleanup
    chainValues = np.unique(slicedCoords.chain)
    count=0
    k=0
    dropList=[]
    for i in tqdm(range(len(chainValues))):
        chainValue = chainValues[i]
        resID = slicedCoords[(slicedCoords.chain == chainValue)].resids
        resNum = len(np.unique(resID))
        if resNum < originalLengths[i]:
            dropIndices = slicedCoords[(slicedCoords.chain == chainValue)].index
            slicedCoords.drop(dropIndices, inplace=True)
            dropList.append(i)
            count += 1
        k+=1
        if k == len(origChains):
            k=0
    i=0
    for x in dropList:
        del originalLengths[x-i]
        i+=1
    print("Number of fragmented chains dropped: {}".format(count))

    # Broken residue edge check
    if brokenResCheck:
        edgeCoords = slicedCoords[slicedCoords.z < 10.0]
        if len(edgeCoords) > 0:
            chainAndResidValues = np.vstack((edgeCoords.chain, edgeCoords.resids)).T
            checkValues = np.vstack(tuple(set(map(tuple,chainAndResidValues))))
            count=0
            for i in tqdm(range(len(checkValues))):
                chainValue, residValue = checkValues[i]
                resname = slicedCoords[(slicedCoords.chain == chainValue) & (slicedCoords.resids == residValue)].resname
                lenRes = len(resname)
                valRes = np.unique(resname)
                if valRes != 0:
                    lenCheck = aminoAcidChecks[str(valRes[0]).strip()]
                    if lenRes != lenCheck:
                        dropIndices = slicedCoords[(slicedCoords.chain == chainValue) & (slicedCoords.resids == residValue)].index
                        slicedCoords.drop(dropIndices, inplace=True)
                        count += 1
            print("Number of fragmented residues dropped: {}".format(count))
            # Drop newly fragmented chains
            chainValues = np.unique(slicedCoords.chain)
            dropList=[]
            count=0
            k=0
            for i in tqdm(range(len(chainValues))):
                chainValue = chainValues[i]
                resID = slicedCoords[(slicedCoords.chain == chainValue)].resids
                resNum = len(np.unique(resID))
                if resNum < originalLengths[i]:
                    dropIndices = slicedCoords[(slicedCoords.chain == chainValue)].index
                    slicedCoords.drop(dropIndices, inplace=True)
                    dropList.append(i)
                    count += 1
                k+=1
                if k == len(origChains):
                    k=0
            print("Number of fragmented chains dropped: {}".format(count))
            i=0
            for x in dropList:
                del originalLengths[x-i]
                i+=1
        else:
            print("No fragmented residues. Nice!")

    # Rename chains
    chainValues = np.unique(slicedCoords.chain)
    chainDictionary = createDictionary.createChainDictionary(chainValues)
    slicedCoords.chain = slicedCoords.chain.replace(chainDictionary)
    print(np.unique(slicedCoords.chain))
    if len(np.unique(slicedCoords.chain)) > 60:
        print("Heads up! There are more than 60 unique chains. Some IDs may not be recognized.")

    # Adjust z for buffer zones
    slicedCoords["z"] = slicedCoords.z.astype(float) - np.min(slicedCoords.z.astype(float))

    # Update index
    slicedCoords["idx"] = range(1,len(slicedCoords)+1)

    # Final formatting and save
    unformattedCoords = copy.deepcopy(slicedCoords)
    slicedCoords = formatPDB.formatPDB(slicedCoords)
    try:
        savetxt(f'output/{output_title}-sliced-w{weight}.pdb', slicedCoords, fmt='%s', delimiter='')
    except:
        print("Couldn't save the sliced capsid! Likely too many unique chains.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=str, default="9jjh", help="PDB filename (without input/)")
    parser.add_argument("--symmetry", type=int, default=5, help="Symmetry of virus")
    parser.add_argument("--plot", action='store_true', help="Enable plotting")
    parser.add_argument("--refindex", type=int, nargs="+", default=[1058], help="Reference indices for plane")
    parser.add_argument("--weight", type=float, default=0.5, help="Slicing weight from 0 to 1")
    parser.add_argument("--auto", action="store_true", help="Automatically detect symmetry axis and reference atom")
    parser.add_argument("--axis-index", type=int, default=0, help="Which candidate axis to use (0=best; ROI-ranked when --roi-selection is set)")
    parser.add_argument("--roi-selection", type=str, default=None, help="MDAnalysis selection for ROI-aware automatic axis selection")
    parser.add_argument("--roi-frame", type=int, default=0, help="Frame/MODEL containing the ROI copy described by --roi-selection")
    parser.add_argument("--list-axes", action="store_true", help="Print candidate axes and exit unless --auto is also set")
    args = parser.parse_args()
    main(args)
