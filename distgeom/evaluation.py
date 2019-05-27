import sys

import numpy as np

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, TorsionFingerprints

isDebug = False

refFileName = '../data/platinum_dataset_2017_01.sdf'
predFileName = sys.argv[1]

# Calculate RMSD with RDKit
refSpl = Chem.SDMolSupplier(refFileName)
predSpl = Chem.SDMolSupplier(predFileName)

cnt = 0
rmsd_list = []
same = 0
for ref, pred in zip(refSpl, predSpl):
    refEntry = ref.GetProp('_Name')
    if pred is None:  # in case of failure
        rmsd = ''
    else:
        try:
            rmsd = AllChem.GetBestRMS(ref, pred)
            rmsd_list.append(rmsd)
        except:
            rmsd = ''
    refInchikey = Chem.inchi.MolToInchiKey(ref)
    predInchikey = Chem.inchi.MolToInchiKey(pred)
    if refInchikey == predInchikey:
        isSame = 'T'
        same += 1
    else:
        isSame = 'F'
    print("{},{},{}".format(refEntry, isSame, rmsd))
    cnt += 1
    if cnt > 3000:
        break
print("same: {}/{}, rmsd: {}".format(same, cnt, np.mean(rmsd_list)))
