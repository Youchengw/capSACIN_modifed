import numpy as np
from numpy import savetxt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import MDAnalysis as md
from tqdm import tqdm
import pandas as pd
import definePlane
import formatPDB
import createDictionary
import copy

def split_dataframe(df, chunk_size): 
    chunks = list()
    num_chunks = len(df) // chunk_size + 1
    for i in range(num_chunks):
        chunks.append(df[i*chunk_size:(i+1)*chunk_size])
    return chunks

#%% USER INPUT!!

output_title = "9jjh"
weight=0.5 # how strict should cut-off be? 0 = no slice, 0.5 = half capsid, 1 = no output

brokenResCheck = True
rotationCheck = True
Plot = False

# indicesVMD = [153841,205513,149535] # Indices (from VMD) for plane generation
# indicesVMD = [138748,147360,134442] # PPV 5-fold
# indicesVMD = [153841,205513,149535] # PPV 3-fold
# indicesVMD = [204523,157603,148552] # PPV 2-fold
# indicesVMD = [5315,14028,18378] # CPV 5-fold
# indicesVMD = [44402] # 3j31
# indicesVMD = [47090] # 2gsy
# indicesVMD = [31968, 31975, 23773] # 3jb8
# indicesVMD = [3903] # 3jb8
# indicesVMD = [5686] # 8hbg
# indicesVMD = [1495] # 3j6r
# indicesVMD = [6020] # 7qwz
# indicesVMD = [1558] # 8hbj

# indicesVMD=[4073] # ppv 2-fold
# indicesVMD=[3131] # ppv 3-fold
# indicesVMD=[959] # ppv 5-fold

# indicesVMD = [936] # 1dzl
# indicesVMD = [821] # 1wcd
# indicesVMD = [996] # 2buk
# indicesVMD = [1210] # 2ztn
# indicesVMD = [601] # 3r0r
# indicesVMD = [928] # 3ra2
# indicesVMD = [591] # 4oq8
# indicesVMD = [592] # 5cw0
# indicesVMD = [1585] # 6jja
# indicesVMD = [2314] # 8des
# indicesVMD = [99] # 9clj
indicesVMD = [1058] # 9jjh
symmetry = 5

gro=f"input/{output_title}.pdb"
traj=f"input/{output_title}.pdb"
u=md.Universe(gro, traj)

origChains = np.unique(u.select_atoms("protein").chainIDs)
nMonomers = len(origChains)

#Obtain reference position
refPos = u.select_atoms(f"protein and index {indicesVMD[0]}")
pointA = refPos.positions[0]

#Create dictionary -- important for removing broken residues
aminoAcidChecks = createDictionary.createDictionary(u)

#Build Virus
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

#Get PDB information
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

#%% Stack the coords
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

#%% Find symmetric reference points
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
 
#%% Alignment infromation 
# indices = np.asarray(indicesVMD) + 1
# pointA, pointB, pointC = chains[chains.idx == indices[0]], chains[chains.idx == indices[1]], chains[chains.idx == indices[2]]
# points=pd.concat((pointA, pointB, pointC))
points = pd.DataFrame(np.vstack((pointA, pointB, pointC)), columns=['x', 'y', 'z'])
x,y,z = points.x, points.y, points.z
selectPoints = np.vstack((x,y,z)).T
com = (np.mean(x), np.mean(y), np.mean(z))
coords = coords - com # center to user-inputed points
selectPoints = selectPoints - com
x, y, z = selectPoints[:,0], selectPoints[:,1], selectPoints[:,2]

#%% Alignment procedure
normalVector = definePlane.definePlane(x,y,z)
zAxis = [0,0,1]

#%%
def rotationMatrix(vec1, vec2):
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.clip(np.dot(a, b), -1.0, 1.0)
    s = np.linalg.norm(v)
    if s < 1e-12:
        if c > 0:
            return np.eye(3)
        helper = np.array([1.0, 0.0, 0.0])
        if abs(a[0]) > 0.9:
            helper = np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, helper)
        axis = axis / np.linalg.norm(axis)
        kmat = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        return np.eye(3) + 2.0 * kmat.dot(kmat)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

if rotationCheck:
    A = rotationMatrix(normalVector, zAxis).T
    newCoords = np.asarray((coords @ A),dtype=float)
    newSelectPoints = np.asarray((selectPoints @ A),dtype=float)
else:
    newCoords = coords
    newSelectPoints = selectPoints

#%%Plot
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

#%%Plot
com = (np.mean(newCoords[:,0]), np.mean(newCoords[:,1]), np.mean(newCoords[:,2]))
centerDist = np.sqrt((newCoords[:,0]-com[0])**2 + (newCoords[:,1]-com[1])**2 + (newCoords[:,2]-com[2])**2)
       
if Plot: 
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(newCoords[0:-1:50,0], newCoords[0:-1:50,1], newCoords[0:-1:50,2], c=centerDist[0:-1:50])
    ax.scatter(newSelectPoints[:,0], newSelectPoints[:,1], newSelectPoints[:,2], color="red", s=90)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.tight_layout()
    plt.show()

#%% Center to COM of new coordinates
newCoords = newCoords - com
newCoords = newCoords - [0,0,np.min(newCoords[:,2])] # Move minimum z to 0

# added 01/05/2024; need to move virus so x and y do not cross -100
newCoords[:,0] = newCoords[:,0] - np.min(newCoords[:,0])
newCoords[:,1] = newCoords[:,1] - np.min(newCoords[:,1])
# z will still be cut-off for now... but I think that's okay. This is just a check anyway.

if Plot:
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(newCoords[0:-1:50,0], newCoords[0:-1:50,1], newCoords[0:-1:50,2], c=centerDist[0:-1:50])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.tight_layout()
    plt.show()

#%% save after alignment, prior to slice
newCoords = pd.DataFrame(newCoords, columns=["x","y","z"])

# update coordinates in created .pdb
chains['x'] = newCoords['x'].apply('{0:.3f}'.format)
chains['y'] = newCoords['y'].apply('{0:.3f}'.format)
chains['z'] = newCoords['z'].apply('{0:.3f}'.format)
chains['b'] = chains['b'].astype('float').apply('{0:.2f}'.format)
chains['occ'] = chains['occ'].astype('float').apply('{0:.2f}'.format)

chains = chains.astype('string')

# update to generate proper .pdb formatting
# ljust --> add blank spaces to right
# rjust --> add blank spaces to left
chains = formatPDB.formatPDB(chains)

# try:
#     savetxt('{}/{}.pdb'.format(opath, output_title), chains, fmt='%s', delimiter='')
# except:
#     print("Couldn't save the full capsid!")
#     print("Likely, there are too many unique chain IDs.")
#     print("This will not cause an issue with slicing.")

chains["x"] = chains["x"].astype('float')
chains["y"] = chains["y"].astype('float')
chains["z"] = chains["z"].astype('float')

#%% time to slice!
cutoff = np.max(chains.z)*weight
originalLen = len(chains)

slicedChains = np.unique(chains.chain[chains["z"] >= cutoff])

originalLengths=[]
for x in slicedChains:
    originalLengths.append(len(np.unique(chains.resids[chains.chain == x])))
slicedCoords = chains[chains["z"] >= cutoff]

slicedCoords.z = slicedCoords.z-(np.min(slicedCoords.z)) # Move sliced capsid such that min z is 0

centerDist = np.sqrt((slicedCoords.x-com[0])**2 + (slicedCoords.y-com[1])**2 + (slicedCoords.z-com[2])**2)

if Plot:
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(slicedCoords.x, slicedCoords.y, slicedCoords.z, c=centerDist)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_zlim(-120,120)
    plt.tight_layout()
    plt.show()

#%% Remove broken chains
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

#%% remove residues that were broken at the boundary
if brokenResCheck == True:
    # Then drop broken residues
    edgeCoords = slicedCoords[slicedCoords.z < 10.0] # Grab positions near slice edge
    checkEdge = len(edgeCoords)
    if checkEdge > 0:
        chainValues = np.unique(edgeCoords.chain)
        residValues = np.unique(edgeCoords.resids)
        chainAndResidValues = np.vstack((edgeCoords.chain, edgeCoords.resids)).T
        checkValues = np.vstack(tuple(set(map(tuple,chainAndResidValues)))) #obtain all unique combinations of chain and resid left in truncated capsid
        
        count=0
        for i in tqdm(range(len(checkValues))):
            chainValue = checkValues[i][0]
            residValue = checkValues[i][1]
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
        
        ### Then drop newly fragmented chains
        dropList=[]
        chainValues = np.unique(slicedCoords.chain)
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

#%% rename chains to be sensible
chainValues = np.unique(slicedCoords.chain)
chainDictionary = createDictionary.createChainDictionary(chainValues)
slicedCoords.chain = slicedCoords.chain.replace(chainDictionary)
print(np.unique(slicedCoords.chain))
if len(np.unique(slicedCoords.chain)) > 60:
    print("Heads up! There are more than 60 unique chains. There may be some IDs that your software of choice can't recognize.")
# slicedCoords = slicedCoords[slicedCoords.chain != "B"]

#%% update for proper .pdb formatting

#one last translation to make definition of buffer zones straightforward
#will always be 0 < z < 5 for fully constrained region
#and will be 5 < z < 15 for C-alpha restrained region
slicedCoords["z"] = slicedCoords.z.astype(float) - np.min(slicedCoords.z.astype(float))

#update index
idList=[]
for i in range(len(slicedCoords)):
    idList.append(i+1)
slicedCoords["idx"] = idList
unformattedCoords = copy.deepcopy(slicedCoords)
slicedCoords = formatPDB.formatPDB(slicedCoords)

# save after slice
try:
    savetxt('output/{}-sliced-w{}.pdb'.format(output_title, weight), slicedCoords, fmt='%s', delimiter='')
except:
    print("Couldn't save the sliced capsid!")
    print("Likely, there are too many unique chain IDs.")
    print("Try adjusting the slicing threshold.")
