from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles("O=C([O-])[C@H]1[C@@H](C=CC=C1C(=O)CCC(=O)[O-])O")
#mol = Chem.AddHs(mol)
#AllChem.EmbedMolecule(mol)
suppl = Chem.SDMolSupplier('164_4MYS_A.sdf')
mol = suppl[0]
print(Chem.MolToSmiles(mol))
print(Chem.inchi.MolToInchiKey(mol))
for a in mol.GetAtoms():
    if a.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
        print("center: ", a.GetIdx())
        print("neighbors: ", ["{}->{}".format(n.GetIdx(), n.GetAtomicNum()) for n in a.GetNeighbors()])
        print("CCW")
    if a.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
        print("center: ", a.GetIdx())
        print("neighbors: ", ["{}->{}".format(n.GetIdx(), n.GetAtomicNum()) for n in a.GetNeighbors()])
        print("CW")
