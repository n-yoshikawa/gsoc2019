import sys
import numpy as np
import pickle
import rdkit.DistanceGeometry as DG
import rdkit
from rdkit import Chem
from rdkit.Chem import rdDistGeom, ChemicalForceFields, rdMolAlign

print(rdkit.__version__)

with open('fragments.pickle', mode='rb') as f:
    db = pickle.load(f)

mol = Chem.MolFromSmiles("CC(C)Cn1c(=O)n(C)c(=O)c2c(-c3ccncc3)n(Cc3cccc4ccccc34)nc21")
canonical_order = list(Chem.rdmolfiles.CanonicalRankAtoms(mol, breakTies=True))
print(canonical_order)

bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)

# Cut input molecule by rotatable bonds
RotatableBond = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
print(Chem.MolToSmiles(mol))
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
for fragment in fragments:
    if fragment.GetNumHeavyAtoms() < 5:
        continue
    fragment_smiles = Chem.MolToSmiles(fragment)
    if fragment_smiles not in db:
        continue
    print(fragment_smiles)
    distMat = db[fragment_smiles]
    for match in mol.GetSubstructMatches(fragment):
        print("match:", match)
        print("canonical:", [canonical_order[i] for i in match])
        for i, p in enumerate(match):
            for j, q in enumerate(match):
                if i == j:
                    continue
                elif i < j:
                    bm[p, q] = distMat[i, j] + 0.0001
                else:
                    bm[p, q] = distMat[i, j] - 0.0001

DG.DoTriangleSmoothing(bm)
ps = rdDistGeom.EmbedParameters()
ps.useRandomCoords = True
ps.SetBoundsMat(bm)
ps.randomSeed = 0xf00d
rdDistGeom.EmbedMolecule(mol, ps)
w = Chem.SDWriter('result.sdf')
w.write(mol)
