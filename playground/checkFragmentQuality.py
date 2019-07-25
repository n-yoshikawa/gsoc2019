import copy
import os
import sys
import numpy as np
np.set_printoptions(precision=4)
import pickle
import rdkit.DistanceGeometry as DG
import rdkit
from rdkit import Chem
from rdkit.Chem import rdDistGeom, ChemicalForceFields, rdMolAlign
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt

print(rdkit.__version__)

with open('fragments-withHs.pickle', mode='rb') as f:
    db = pickle.load(f)

suppl = Chem.SDMolSupplier('platinum_dataset_2017_01.sdf')
cnt_fragment = 0
cnt_found = 0
cnt_rmsd_error = 0
rmsd_list = []
for mol in suppl:
    smiles = Chem.MolToSmiles(mol)
    mol = Chem.AddHs(mol)
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

    for fragment in fragments:
        # Skip small fragments
        if fragment.GetNumHeavyAtoms() < 5:
            continue
        cnt_fragment += 1
        inchikey = Chem.inchi.MolToInchiKey(fragment)
        if os.path.exists('{}.sdf'.format(inchikey)):
            cnt_found += 1
            db_sdf = Chem.SDMolSupplier('{}.sdf'.format(inchikey))
            fragment_from_db = db_sdf[0]
            fragment = Chem.AddHs(fragment)
            fragment_from_db = Chem.AddHs(fragment_from_db)
            try:
                rmsd = AllChem.GetBestRMS(fragment, fragment_from_db)
                rmsd_list.append(rmsd)
                print('{},{}'.format(inchikey, rmsd))
            except:
                cnt_rmsd_error += 1
                continue
print(cnt_fragment, cnt_found, cnt_rmsd_error)
plt.hist(rmsd_list, bins=16)
plt.show()
