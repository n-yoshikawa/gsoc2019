import io
import sys
import pickle
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import rdDistGeom
np.set_printoptions(precision=4)
print(rdkit.__version__)

suppl = Chem.SDMolSupplier(sys.argv[1])
fragment_list = {}
for mol in suppl:
    if mol is None:
        continue
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
    try:
        fragments = Chem.rdmolops.GetMolFrags(rwmol.GetMol(), asMols=True)
    except:
        continue
    # Generate distance matrix from fragments
    for fragment in fragments:
        # Skip small fragments
        if fragment.GetNumHeavyAtoms() < 5:
            continue
        smiles = Chem.MolToSmiles(fragment)
        print(smiles)
        atomOrder = fragment.GetPropsAsDict(includePrivate=True,includeComputed=True)['_smilesAtomOutputOrder']
        
        fragment = Chem.RenumberAtoms(fragment,atomOrder)
        
        conf = fragment.GetConformer()
        N = fragment.GetNumAtoms()
        distMat = np.zeros((N, N))
        for atom1 in fragment.GetAtoms():
            for atom2 in fragment.GetAtoms():
                i = atom1.GetIdx()
                j = atom2.GetIdx()
                distMat[i, j] = (conf.GetAtomPosition(i)-conf.GetAtomPosition(j)).Length()
        print(smiles)
        fragment_list[smiles] = distMat
# Save distance matrix by pickle
with open('fragments.pickle', mode='wb') as f:
    pickle.dump(fragment_list, f)
