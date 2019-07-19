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

w = Chem.SDWriter('result.sdf')
with open(sys.argv[1], "r") as f:
    for line in f:
        smiles, entry = line.split()
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp("_Name", entry)
        print("Processing:", smiles)
        
        mol = Chem.AddHs(mol)
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
            #print(fragment_smiles)
            cfragment = Chem.MolFromSmiles(fragment_smiles)
            distMat = db[fragment_smiles]
            #print("non canonical:", mol.GetSubstructMatches(fragment))
            #print("canonical:", mol.GetSubstructMatches(cfragment))
            for match in mol.GetSubstructMatches(cfragment):
                containProcessed = False
                for a in match:
                    if a in processed:
                        containProcessed = True
                        break
                if containProcessed:
                    #print("duplicated!")
                    continue
                processed.extend(match)
                for i, p in enumerate(match):
                    for j, q in enumerate(match):
                        if i == j:
                            continue
                        elif p < q:
                            #print("Update {}, {}".format(p, q))
                            bm[p, q] = distMat[i, j] + 0.01
                        else:
                            bm[p, q] = distMat[i, j] - 0.01
        
        #print("Pre-smoothing check")
        for i in range(len(bm)):
            for j in range(len(bm)):
                if i < j:
                    if bm[j, i] > bm[j, i]:
                        print("({}, {}): {} < x < {}".format(i, j, bm[j, i], bm[i, j]))
                        print("({}, {}): {} < x < {} (original)".format(i, j, bm_org[j, i], bm_org[i, j]))
                else:
                    if bm[i, j] > bm[j, i]:
                        print("({}, {}): {} < x < {}".format(i, j, bm[i, j], bm[j, i]))
                        print("({}, {}): {} < x < {} (original)".format(i, j, bm_org[i, j], bm_org[j, i]))
        #DG.DoTriangleSmoothing(bm)
        #print("Post-smoothing check")
        #for i in range(len(bm)):
        #    for j in range(len(bm)):
        #        if i < j:
        #            if bm[j, i] > bm[i, j]:
        #                print("({}, {}): {} < x < {}".format(i, j, bm[j, i], bm[i, j]))
        #                print("({}, {}): {} < x < {} (original)".format(i, j, bm_org[j, i], bm_org[i, j]))
        #        else:
        #            if bm[i, j] > bm[j, i]:
        #                print("({}, {}): {} < x < {}".format(i, j, bm[i, j], bm[j, i]))
        #                print("({}, {}): {} < x < {} (original)".format(i, j, bm_org[i, j], bm_org[j, i]))
        ps = rdDistGeom.ETKDG()
        ps.useRandomCoords = True
        ps.SetBoundsMat(bm)
        ps.randomSeed = 0xf00d
        try:
            rdDistGeom.EmbedMolecule(mol, ps)
        except:
            print("Failed to generate coordinates")
            sys.exit()
        ps = rdDistGeom.ETKDG()
        ps.useRandomCoords = True
        ps.randomSeed = 0xf00d
        rdDistGeom.EmbedMolecule(mol, ps)
        w.write(mol)
