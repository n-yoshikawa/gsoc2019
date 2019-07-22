import copy
import sys
import numpy as np
np.set_printoptions(precision=4)
import pickle
import rdkit.DistanceGeometry as DG
import rdkit
from rdkit import Chem
from rdkit.Chem import rdDistGeom, ChemicalForceFields, rdMolAlign

print(rdkit.__version__)

with open('fragments.pickle', mode='rb') as f:
    db = pickle.load(f)

smiles = "s1c(c(c(c1c1cc(ccc1)NC1CCCCC1)Br)OCC(=O)[O-])C(=O)[O-]"
mol = Chem.MolFromSmiles(smiles)

bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
bm_org = copy.deepcopy(bm)

# Cut input molecule by rotatable bonds
RotatableBond = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
rwmol = Chem.RWMol(mol)
for begin, end in mol.GetSubstructMatches(RotatableBond):
    rwmol.RemoveBond(begin, end)
    beginAtom = rwmol.GetAtomWithIdx(begin)
    endAtom = rwmol.GetAtomWithIdx(end)
    if beginAtom.GetAtomicNum() != 6 and beginAtom.GetIsAromatic():
        beginAtom.SetNumExplicitHs(1)
        beginAtom.SetNoImplicit(True)
    if endAtom.GetAtomicNum() != 6 and endAtom.GetIsAromatic():
        endAtom.SetNumExplicitHs(1)
        endAtom.SetNoImplicit(True)
fragments = Chem.rdmolops.GetMolFrags(rwmol.GetMol(), asMols=True)

# Set boundary from fragments
processed = []
for fragment in fragments:
    if fragment.GetNumHeavyAtoms() < 5:
        continue
    fragment_smiles = Chem.MolToSmiles(fragment)
    if fragment_smiles not in db:
        continue
    print(fragment_smiles)
    cfragment = Chem.MolFromSmiles(fragment_smiles)
    distMat = db[fragment_smiles]
    print("match:", mol.GetSubstructMatches(cfragment))
    for match in mol.GetSubstructMatches(cfragment):
        containProcessed = False
        for a in match:
            if a in processed:
                containProcessed = True
                break
        if containProcessed:
            continue
        processed.extend(match)
        for i, p in enumerate(match):
            for j, q in enumerate(match):
                if i == j:
                    continue
                elif p < q:
                    bm[p, q] = distMat[i, j] + 0.01
                else:
                    bm[p, q] = distMat[i, j] - 0.01
def ub(i, j):
    if i < j:
        return bm[i, j]
    else:
        return bm[j, i]
def lb(i, j):
    if i < j:
        return bm[j, i]
    else:
        return bm[i, j]

# assure bounds are nonerroneous 
for i in range(len(bm)):
    for j in range(len(bm)):
        assert lb(i, j) <= ub(i, j)
        assert lb(i, j) >= 0.0
DG.DoTriangleSmoothing(bm)
for i in range(len(bm)):
    for j in range(len(bm)):
        assert lb(i, j) <= ub(i, j)
        assert lb(i, j) >= 0.0
