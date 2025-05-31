from ase.build import bulk, make_supercell
from ase.io import write
from icet import ClusterSpace
from icet.tools.structure_generation import generate_sqs
import os
import numpy as np
import time

# Step 1: Define elements
elements = ["Zr", "Nb", "Ti", "V", "Hf"]

# Step 2: Create prototype structure (BCC with supercell)
prototype = bulk("V", crystalstructure="bcc", a=3.3)
P = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])  # 2x2x2 supercell
prototype = make_supercell(prototype, P)

# Step 3: Define cutoff distances
cutoffs = [7.0, 5.0, 4.0]

# Step 4: Define chemical_symbols for each atomic site
chemical_symbols = [list(elements) for _ in range(len(prototype))]

# Debug check
print("len(prototype) =", len(prototype))
print("chemical_symbols =", chemical_symbols)
assert len(chemical_symbols) == len(prototype)

# Step 5: Create ClusterSpace
cluster_space = ClusterSpace(
    structure=prototype,
    chemical_symbols=chemical_symbols,
    cutoffs=cutoffs
)
print("✅ Cluster space created with", len(cluster_space), "clusters")

# Step 6: Define concentrations
target_concentrations = {el: 1 / len(elements) for el in elements}

# Step 7: Create output directory
output_dir = "outputs"
os.makedirs(output_dir, exist_ok=True)
print(f"📁 Output directory ready: {output_dir}")

# Step 8: Generate SQS structures
for seed in [42, 43, 44]:
    try:
        print(f"\n🔁 Generating SQS for seed {seed}...")
        start_time = time.time()
        structure = generate_sqs(
            cluster_space=cluster_space,
            target_concentrations=target_concentrations,
            random_seed=seed,
            n_steps=10000,
            max_size=64
        )
        duration = time.time() - start_time
        print(f"✅ SQS generation complete in {duration:.2f} seconds")
        print(f"🧱 Final structure has {len(structure)} atoms")

        # Save the structure
        output_path = os.path.join(output_dir, f"sqs_seed_{seed}.xyz")
        write(output_path, structure)
        print(f"📦 Saved structure to {output_path}")

    except Exception as e:
        print(f"❌ Error generating seed {seed}: {e}")
