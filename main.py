import numpy as np
from Bio import PDB
import matplotlib.pyplot as plt

pdb_id = '4YWO'
parser = PDB.PDBParser()
structure = parser.get_structure(pdb_id, f"{pdb_id}.pdb")

ca_atoms = []
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.has_id('CA'):
                ca_atoms.append(residue['CA'])

num_ca = len(ca_atoms)
contact_map = np.zeros((num_ca, num_ca))

for i in range(num_ca):
    for j in range(i + 1, num_ca):
        distance = ca_atoms[i] - ca_atoms[j]
        if distance < 8.0:  # prog 8A
            contact_map[i][j] = 1
            contact_map[j][i] = 1

plt.imshow(contact_map, cmap='binary')
plt.title('Mapa kontaktów')
plt.xlabel('Atom CA')
plt.ylabel('Atom CA')
plt.show()
