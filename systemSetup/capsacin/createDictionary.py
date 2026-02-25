import MDAnalysis as md
import itertools

def createDictionary(u):
    aminoAcids = ["ALA","ARG","ASN","ASP","CYS","GLU","GLN",
                  "GLY","HIS","ILE","LEU","LYS","MET","PHE",
                  "PRO","SER","THR","TRP","TYR","VAL"]
    
    counts = []
    for i in aminoAcids:
        try:
            atoms = len(u.select_atoms("protein and resname {} and chainid A".format(i)).split("residue")[0].positions)
        except:
            atoms = 0
        counts.append(atoms)
    
    res = {}
    for key in aminoAcids:
        for value in counts:
            res[key] = value
            counts.remove(value)
            break
    
    return res

def createChainDictionary(chainValues):
    # chainID = [chr(65 + i) for i in range(len(chainValues))]
    # Define character sets in order of priority
    capital_letters = [chr(i) for i in range(65, 91)]  # A-Z
    # lowercase_letters = [chr(i) for i in range(97, 123)]  # a-z
    numbers = [str(i) for i in range(10)]  # 0-9
    symbols = list("!@#$%^&-_=+;:'\",.<>?/\\|`~")
    
    # Chain them together
    all_tokens = list(itertools.chain(capital_letters, numbers, symbols))
    chainID = all_tokens[:len(chainValues)]
    
    res = {}
    for key in chainValues:
        for value in chainID:
            res[key] = value
            chainID.remove(value)
            break
    
    return res
