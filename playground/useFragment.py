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
        print("Processing:", entry, smiles)
        
        
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

        mol = Chem.AddHs(mol)
        bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
        bm_org = copy.deepcopy(bm)
        
        # Set boundary from fragments
        for fragment in fragments:
            if fragment.GetNumHeavyAtoms() < 5:
                continue
            fragment_smiles = Chem.MolToSmiles(fragment)
            cfragment = Chem.MolFromSmiles(fragment_smiles)
            csmiles = Chem.MolToSmiles(cfragment)
            if fragment_smiles not in db:
                continue
            cfragment = Chem.MolFromSmarts(fragment_smiles)
            distMat = db[fragment_smiles]
            #print("non canonical:", mol.GetSubstructMatches(fragment))
            #print("{}, match: {}".format(csmiles.rstrip(), mol.GetSubstructMatches(cfragment)))
            processed = []
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
                            bm[p, q] = distMat[i, j] + 0.1
                            bm[q, p] = max(0, distMat[i, j] - 0.1)

                        else:
                            bm[p, q] = max(0, distMat[i, j] - 0.1)
                            bm[q, p] = distMat[i, j] + 0.1
                    for q in range(len(bm)):
                        if q in match:
                            continue
                        if mol.GetAtomWithIdx(q).GetSymbol() != 'H':
                            if p < q:
                                bm[p, q] = 1000
                                bm[q, p] = 0
                            else:
                                bm[p, q] = 0
                                bm[q, p] = 1000
        #print("original:")
        #for i in range(len(bm)):
        #    for j in range(i+1, len(bm)):
        #        print("({}, {}) {} < x < {}".format(i, j, bm_org[j, i], bm_org[i, j]))
        for i in range(len(bm)):
            for j in range(i+1, len(bm)):
                #print("({}, {}) {} < x < {}".format(i, j, bm[j, i], bm[i, j]))
                if bm[i, j] < bm[j, i]:
                    print("   ***** Assertion failed !! (Before smoothing) *****")
                assert(bm[i, j] >= bm[j, i])

        if DG.DoTriangleSmoothing(bm) == False:
            print("Smoothing failed")

        for i in range(len(bm)):
            for j in range(i+1, len(bm)):
                #print("({}, {}) {} < x < {}".format(i, j, bm[j, i], bm[i, j]))
                if bm[i, j] < bm[j, i]:
                    print("   ***** Assertion failed !! (After smoothing) *****")
                assert(bm[i, j] >= bm[j, i])
        ps = rdDistGeom.EmbedParameters()
        ps.useRandomCoords = True
        ps.SetBoundsMat(bm)
        ps.randomSeed = 0xf00d
        try:
            rdDistGeom.EmbedMolecule(mol, ps)
        except:
            print("Failed to generate coordinates")
            ps = rdDistGeom.ETKDG()
            ps.useRandomCoords = True
            ps.randomSeed = 0xf00d
            rdDistGeom.EmbedMolecule(mol, ps)
        w.write(mol)
