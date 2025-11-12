from rdkit import Chem
from rdkit.Chem import AllChem

smiles = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O"
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

# Try multiple embedding attempts with ETKDGv3
embed_status = AllChem.EmbedMultipleConfs(
    mol,
    numConfs=10,
    params=AllChem.ETKDGv3()
)

if len(embed_status) == 0:
    raise ValueError("Embedding failed for all attempts. Try increasing numConfs or using randomSeed.")

# Optimize all conformers
for conf_id in embed_status:
    AllChem.UFFOptimizeMolecule(mol, confId=conf_id)

# Write all conformers to SDF
w = Chem.SDWriter("estradiol.sdf")
for conf_id in embed_status:
    w.write(mol, confId=conf_id)
w.close()
print(f"SDF file generated with {len(embed_status)} conformers.")

