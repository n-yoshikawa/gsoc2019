import io
import sys
import pickle
import numpy as np
import rdkit
from rdkit import Chem

print(rdkit.__version__)

suppl = Chem.SDMolSupplier(sys.argv[1])
fragment_list = {}
for mol in suppl:
    print(Chem.MolToSmiles(mol))
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
    # Generate distance matrix from fragments
    for fragment in fragments:
        # Skip small fragments
        if fragment.GetNumHeavyAtoms() < 5:
            continue
        conf = fragment.GetConformer()
        N = fragment.GetNumAtoms()
        distMat = np.zeros((N, N))
        for atom1 in fragment.GetAtoms():
            for atom2 in fragment.GetAtoms():
                i = atom1.GetIdx()
                j = atom2.GetIdx()
                distMat[i, j] = (conf.GetAtomPosition(i)-conf.GetAtomPosition(j)).Length()
        smiles = Chem.MolToSmiles(fragment)
        fragment_list[smiles] = distMat
# Save distance matrix by pickle
with open('fragments.pickle', mode='wb') as f:
    pickle.dump(fragment_list, f)

# For debug
# Get canonical order
canonical_order = list(Chem.rdmolfiles.CanonicalRankAtoms(mol, breakTies=True))
print(canonical_order)

# Print canonical order of match
for fragment_smiles in fragment_list:
    print(fragment_smiles)
    fragment = Chem.MolFromSmiles(fragment_smiles)
    for match in mol.GetSubstructMatches(fragment):
        print("match:", match)
        print("canonical:", [canonical_order[i] for i in match])
