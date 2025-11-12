from Bio.PDB import PDBParser
import numpy as np
parser = PDBParser()
structure = parser.get_structure("protein", "4run_chainA_protein.pdb")
coords = []
for atom in structure.get_atoms():
    if atom.get_parent().get_id()[0] == " ":  # ignore heteroatoms
        coords.append(atom.get_coord())
coords = np.array(coords)
centroid = coords.mean(axis=0)
print("Centroid:", centroid)
