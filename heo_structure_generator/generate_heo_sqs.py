import random
import os
from ase.build import bulk
from ase.io import write
from icet import ClusterSpace
from icet.tools.structure_generation import generate_sqs

# Set random seed for reproducibility
random.seed(42)

# Define A-site elements for high-entropy oxide
a_site_elements = ['Mg', 'Ni', 'Co', 'Cu', 'Zn']

# Create output directory
os.makedirs("outputs", exist_ok=True)

# Step 1: Build a rocksalt MgO primitive cell
atoms = bulk('MgO', crystalstructure='rocksalt', a=4.2, cubic=True)

# Step 2: Replace all Mg atoms with a dummy element "X"
for atom in atoms:
    if atom.symbol == 'Mg':
        atom.symbol = 'X'

# Step 3: Set up the ClusterSpace
cluster_space = ClusterSpace(atoms, cutoffs=[5.0], chemical_symbols=a_site_elements)

# Step 4: Define target composition (equal ratio)
target = {element: 1.0 / len(a_site_elements) for element in a_site_elements}

# Step 5: Generate the SQS structure
sqs = generate_sqs(cluster_space, target_occupations={"X": target}, supercell=[2, 2, 2])

# Step 6: Save the generated structure
write("outputs/structure.xyz", sqs)
write("outputs/POSCAR", sqs, format='vasp')

# Print structure info
print(f"Number of atoms: {len(sqs)}")
print(f"Elements in structure: {sqs.get_chemical_symbols()}")
print("✅ Structure saved as 'structure.xyz' and 'POSCAR' in the 'outputs' folder.")
input("Press Enter to exit...")
