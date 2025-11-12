from rdkit import Chem
from rdkit.Chem import AllChem
# Load existing molecule file (testo.mol)
mol = Chem.MolFromMolFile("testo.mol", removeHs=False)
# If hydrogens are missing, add them
mol = Chem.AddHs(mol)
# Generate 3D coordinates 
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.UFFOptimizeMolecule(mol)
# Save as SDF
w = Chem.SDWriter("testosterone.sdf")
w.write(mol)
w.close()
print("Converted testo.mol → testosterone.sdf with hydrogens and 3D coordinates.")
