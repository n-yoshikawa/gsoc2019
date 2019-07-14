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

w = Chem.SDWriter('result.sdf')
with open(sys.argv[1], "r") as f:
    for line in f:
        smiles, entry = line.split()
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp("_Name", entry)
        print(smiles)
        
        bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
        
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
        for fragment in fragments:
            if fragment.GetNumHeavyAtoms() < 5:
                continue
            fragment_smiles = Chem.MolToSmiles(fragment)
            if fragment_smiles not in db:
                continue
            #print(fragment_smiles)
            cfragment = Chem.MolFromSmiles(fragment_smiles)
            distMat = db[fragment_smiles]
            #print("non canonical:", mol.GetSubstructMatches(fragment))
            #print("canonical:", mol.GetSubstructMatches(cfragment))
            for match in mol.GetSubstructMatches(cfragment):
                for i, p in enumerate(match):
                    for j, q in enumerate(match):
                        if i == j:
                            continue
                        elif p < q:
                            bm[p, q] = distMat[i, j] + 0.0001
                        else:
                            bm[p, q] = distMat[i, j] - 0.0001
        
        DG.DoTriangleSmoothing(bm)
        ps = rdDistGeom.EmbedParameters()
        ps.useRandomCoords = True
        ps.SetBoundsMat(bm)
        ps.randomSeed = 0xf00d
        rdDistGeom.EmbedMolecule(mol, ps)
        w.write(mol)
