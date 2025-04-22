import random
from ase.build import bulk
from ase.io import write
import os

# Set random seed
random.seed(42)

# A-site elements for HEO
a_site_elements = ['Mg', 'Ni', 'Co', 'Cu', 'Zn']

# Create output directory
os.makedirs("outputs", exist_ok=True)

# Create MgO rocksalt and expand to 2x2x2
atoms = bulk('MgO', crystalstructure='rocksalt', a=4.2, cubic=True)
atoms *= (2, 2, 2)

# Replace Mg with random A-site elements
for atom in atoms:
    if atom.symbol == 'Mg':
        atom.symbol = random.choice(a_site_elements)

# Print result
print(f"Number of atoms: {len(atoms)}")
print(f"Composition: {atoms.get_chemical_symbols()}")

# Save structure to file
write("outputs/structure.xyz", atoms)
write("outputs/POSCAR", atoms, format='vasp')

print("Structure saved as 'structure.xyz' and 'POSCAR' in the 'outputs' folder.")
input("Press Enter to exit...")
